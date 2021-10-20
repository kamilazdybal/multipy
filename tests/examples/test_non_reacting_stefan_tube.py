import unittest
import numpy as np
import multipy
from scipy.sparse.linalg import expm
from scipy.integrate import odeint
from scipy.optimize import fmin

################################################################################
################################################################################
####
####    Class: NonReactingStefanTube
####
################################################################################
################################################################################

composition = multipy.Composition()
velocity = multipy.Velocity()
transform = multipy.Transform()
check = multipy.Check()

n_species = 3
species_names = ['Acetone', 'Methanol', 'Air']
T = 328.5
p = 101325
mixture_molar_density = composition.mixture_molar_density(T, p)
species_molar_masses = np.array([[58.08/1000], [32.04/1000], [28.9628/1000]])
D12 = 8.48 / 1000**2
D21 = D12
D13 = 13.72 / 1000**2
D31 = D13
D23 = 19.91 / 1000**2
D32 = D23
D_binary = np.array([[0, D12, D13],
                     [D21, 0, D23],
                     [D31, D32, 0]])

length = 0.238

def phi_Phi(N, D_binary, mixture_molar_density, length):

    Phi = np.zeros((n_species-1,n_species-1))
    phi = np.zeros((n_species-1,1))

    for i in range(0,n_species-1):
        for j in range(0,n_species-1):

            if i == j:

                summed_terms = 0
                for k in range(0, n_species):
                    if k != i:
                        summed_terms = summed_terms + N[k] / (mixture_molar_density * D_binary[i,k]/length)

                Phi[i,i] = N[i] / (mixture_molar_density * D_binary[i,n_species-1]/length) + summed_terms

            else:

                Phi[i,j] = N[i] * (1/(mixture_molar_density * D_binary[i,n_species-1]/length) - 1/(mixture_molar_density * D_binary[i,j]/length))

        phi[i] = - N[i] / (mixture_molar_density * D_binary[i,n_species-1]/length)

    return (phi, Phi)

def b_tilde_B_tilde(omega, n, mixture_molar_density, species_molar_masses, binary_diffusion_coefficients, length):

    omega_acetone, omega_methanol = omega
    omega_air = 1 - omega_acetone - omega_methanol
    species_mass_fractions = np.array([[omega_acetone],[omega_methanol],[omega_air]])

    B_tilde = np.zeros((2,2))
    b_tilde = np.zeros((2,1))

    Mt = composition.mixture_molar_mass(species_mass_fractions, 'mass', species_molar_masses)

    for i in range(0,2):
        for j in range(0,2):

            if i == j:

                summed_terms = 0
                for k in range(0, 3):
                    if k != i:
                        summed_terms = summed_terms + n[k] / (species_molar_masses[k] * binary_diffusion_coefficients[i,k]/length)

                B_tilde[i,i] = Mt / (mixture_molar_density * species_molar_masses[i]) * (n[i]/(species_molar_masses[i] * binary_diffusion_coefficients[i,2]/length) + summed_terms)

            else:

                B_tilde[i,j] = (n[i] * Mt / (mixture_molar_density * species_molar_masses[i] * species_molar_masses[j])) * (1/(binary_diffusion_coefficients[i,2]/length) - 1/(binary_diffusion_coefficients[i,j]/length))

        b_tilde[i] = - n[i] / (mixture_molar_density * species_molar_masses[i] * binary_diffusion_coefficients[i,2]/length)

    return (b_tilde, B_tilde)

def phi_tilde_Phi_tilde(b_tilde, B_tilde, J):

    phi = np.dot(J[:,:,0], b_tilde)

    Phi = np.dot(J[:,:,0], B_tilde)

    return(phi, Phi)

def stefan_tube_mass_fractions_numerical(n_points=1000):

    eta_coordinates = np.linspace(0,1,n_points)
    M_at_initial_condition = 0.319 * species_molar_masses[0,0] + 0.528 * species_molar_masses[1,0] + (1 - 0.319 - 0.528) * species_molar_masses[2,0]
    initial_condition = np.array([species_molar_masses[0,0] * 0.319 / M_at_initial_condition, species_molar_masses[1,0] * 0.528 / M_at_initial_condition])
    boundary_condition = np.array([0,0,1])

    def RHS_ODE(omega, time_vector, n, mixture_molar_density, species_molar_masses, binary_diffusion_coefficients, length):

        omega_acetone, omega_methanol = omega
        omega_air = 1 - omega_acetone - omega_methanol
        species_mass_fractions = np.array([[omega_acetone],[omega_methanol],[omega_air]])

        (b_tilde, B_tilde) = b_tilde_B_tilde(omega, n, mixture_molar_density, species_molar_masses, binary_diffusion_coefficients, length)

        J = transform.species_gradients_mole_to_mass(species_mass_fractions, species_molar_masses)

        (phi_tilde, Phi_tilde) = phi_tilde_Phi_tilde(b_tilde, B_tilde, J)

        gradient_omega_acetone = Phi_tilde[0,0] * omega_acetone + Phi_tilde[0,1] * omega_methanol + phi_tilde[0]
        gradient_omega_methanol = Phi_tilde[1,0] * omega_acetone + Phi_tilde[1,1] * omega_methanol + phi_tilde[1]

        dxdeta = np.array([gradient_omega_acetone, gradient_omega_methanol])

        return dxdeta.ravel()

    def error_function(fluxes):

        (n1, n2) = fluxes
        n = np.array([n1, n2, 0])

        numerical_solution = odeint(RHS_ODE, initial_condition, eta_coordinates, args=(n, mixture_molar_density, species_molar_masses, D_binary, length))

        mass_fractions_numerical = np.vstack((numerical_solution[:,0], numerical_solution[:,1], 1 - numerical_solution[:,0] - numerical_solution[:,1]))

        # Calculate an error between the current mass fractions at the boundary and the boundary condition:
        error = np.linalg.norm(boundary_condition - mass_fractions_numerical[:,-1])

        return error

    minimum = fmin(func=error_function, x0=(0.00001, 0.00001), xtol=10**-8, ftol=10**-8, maxiter=10000, disp=False)

    n1 = minimum[0]
    n2 = minimum[1]
    n3 = 0
    n = np.array([n1, n2, n3])

    numerical_solution = odeint(RHS_ODE, initial_condition, eta_coordinates, args=(n, mixture_molar_density, species_molar_masses, D_binary, length))

    species_mass_fractions = np.vstack((numerical_solution[:,0], numerical_solution[:,1], 1 - numerical_solution[:,0] - numerical_solution[:,1]))

    composition.set_species_mass_fractions = species_mass_fractions

    return species_mass_fractions

def stefan_tube_mole_fractions_numerical(n_points=1000):

    eta_coordinates = np.linspace(0,1,n_points)
    initial_condition = np.array([0.319, 0.528])
    boundary_condition = np.array([0,0,1])

    def RHS_ODE(x, time_vector, N, D_binary, mixture_molar_density, length):

        x1, x2 = x

        (phi, Phi) = phi_Phi(N, D_binary, mixture_molar_density, length)

        gradient_x_acetone = Phi[0,0] * x1 + Phi[0,1] * x2 + phi[0]
        gradient_x_methanol = Phi[1,0] * x1 + Phi[1,1] * x2 + phi[1]

        dxdeta = np.array([gradient_x_acetone, gradient_x_methanol])

        return dxdeta.ravel()

    def error_function(fluxes):

        (N1, N2) = fluxes

        N = np.array([N1, N2, 0])

        numerical_solution = odeint(RHS_ODE, initial_condition, eta_coordinates, args=(N, D_binary, mixture_molar_density, length))

        mole_fractions_numerical = np.vstack((numerical_solution[:,0], numerical_solution[:,1], 1 - numerical_solution[:,0] - numerical_solution[:,1]))

        # Calculate an error between the current mole fractions at the boundary and the boundary condition:
        error = np.linalg.norm(boundary_condition - mole_fractions_numerical[:,-1])

        return error

    minimum = fmin(func=error_function, x0=(0.001, 0.001), xtol=10**-8, ftol=10**-8, disp=False)

    N1 = minimum[0]
    N2 = minimum[1]
    N3 = 0

    N = np.array([N1, N2, N3])

    numerical_solution = odeint(RHS_ODE, initial_condition, eta_coordinates, args=(N, D_binary, mixture_molar_density, length))

    species_mole_fractions = np.vstack((numerical_solution[:,0], numerical_solution[:,1], 1 - numerical_solution[:,0] - numerical_solution[:,1]))

    composition.set_species_mole_fractions = species_mole_fractions

    return species_mole_fractions

def stefan_tube_mole_fractions_analytic(n_points=1000):

    eta_coordinates = np.linspace(0,1,n_points)
    initial_condition = np.array([0.319, 0.528])
    boundary_condition = np.array([0,0,1])

    def RHS_ODE(x, time_vector, N, D_binary, mixture_molar_density, length):

        x1, x2 = x

        (phi, Phi) = phi_Phi(N, D_binary, mixture_molar_density, length)

        gradient_x_acetone = Phi[0,0] * x1 + Phi[0,1] * x2 + phi[0]
        gradient_x_methanol = Phi[1,0] * x1 + Phi[1,1] * x2 + phi[1]

        dxdeta = np.array([gradient_x_acetone, gradient_x_methanol])

        return dxdeta.ravel()

    def error_function(fluxes):

        (N1, N2) = fluxes

        N = np.array([N1, N2, 0])

        numerical_solution = odeint(RHS_ODE, initial_condition, eta_coordinates, args=(N, D_binary, mixture_molar_density, length))

        mole_fractions_numerical = np.vstack((numerical_solution[:,0], numerical_solution[:,1], 1 - numerical_solution[:,0] - numerical_solution[:,1]))

        # Calculate an error between the current mole fractions at the boundary and the boundary condition:
        error = np.linalg.norm(boundary_condition - mole_fractions_numerical[:,-1])

        return error

    minimum = fmin(func=error_function, x0=(0.001, 0.001), xtol=10**-8, ftol=10**-8, disp=False)

    N1 = minimum[0]
    N2 = minimum[1]
    N3 = 0

    N = np.array([N1, N2, N3])

    (phi, Phi) = phi_Phi(N, D_binary, mixture_molar_density, length)

    species_mole_fractions = np.zeros((3,n_points))

    for i, eta in enumerate(eta_coordinates):

        species_mole_fractions[0:2,i:i+1] = np.dot(expm(Phi*eta), initial_condition[:,None]) + np.dot((expm(Phi*eta) - np.identity(2)), np.dot(np.linalg.inv(Phi), phi))
        species_mole_fractions[2,i] = 1 - np.sum(species_mole_fractions[0:2,i])

    composition.set_species_mole_fractions = species_mole_fractions

    return species_mole_fractions

class NonReactingStefanTube(unittest.TestCase):

    def test__Mole_fractions_numerical_computation(self):

        for n_points in [10, 100]:

            try:
                species_mole_fractions = stefan_tube_mole_fractions_numerical(n_points=n_points)
                (n_species, n_observations) = np.shape(species_mole_fractions)
                self.assertTrue(n_species==3)
                self.assertTrue(n_observations==n_points)

                acetone_tube_outlet = round(species_mole_fractions[0,-1], 5)
                methanol_tube_outlet = round(species_mole_fractions[1,-1], 5)
                air_tube_outlet = round(species_mole_fractions[2,-1],5)
                self.assertTrue(acetone_tube_outlet==0)
                self.assertTrue(methanol_tube_outlet==0)
                self.assertTrue(air_tube_outlet==1)

                acetone_liquid_surface = round(species_mole_fractions[0,0], 5)
                methanol_liquid_surface = round(species_mole_fractions[1,0], 5)
                self.assertTrue(acetone_liquid_surface==0.319)
                self.assertTrue(methanol_liquid_surface==0.528)

                idx = check.sum_of_species_fractions(species_mole_fractions, tolerance=0.000001, verbose=False)
                (idx_below_zero, idx_above_one) = check.range_of_species_fractions(species_mole_fractions, tolerance=0.000001, verbose=False)
                self.assertTrue(len(idx)==0)
                self.assertTrue(len(idx_below_zero)==0)
                self.assertTrue(len(idx_above_one)==0)

            except Exception:
                self.assertTrue(False)

    def test__Mole_fractions_analytic_computation(self):

        for n_points in [10, 100]:

            try:
                species_mole_fractions = stefan_tube_mole_fractions_analytic(n_points=n_points)
                (n_species, n_observations) = np.shape(species_mole_fractions)
                self.assertTrue(n_species==3)
                self.assertTrue(n_observations==n_points)

                acetone_tube_outlet = round(species_mole_fractions[0,-1], 5)
                methanol_tube_outlet = round(species_mole_fractions[1,-1], 5)
                air_tube_outlet = round(species_mole_fractions[2,-1],5)
                self.assertTrue(acetone_tube_outlet==0)
                self.assertTrue(methanol_tube_outlet==0)
                self.assertTrue(air_tube_outlet==1)

                acetone_liquid_surface = round(species_mole_fractions[0,0], 5)
                methanol_liquid_surface = round(species_mole_fractions[1,0], 5)
                self.assertTrue(acetone_liquid_surface==0.319)
                self.assertTrue(methanol_liquid_surface==0.528)

                idx = check.sum_of_species_fractions(species_mole_fractions, tolerance=0.000001, verbose=False)
                (idx_below_zero, idx_above_one) = check.range_of_species_fractions(species_mole_fractions, tolerance=0.000001, verbose=False)
                self.assertTrue(len(idx)==0)
                self.assertTrue(len(idx_below_zero)==0)
                self.assertTrue(len(idx_above_one)==0)

            except Exception:
                self.assertTrue(False)

    def test__Mass_fractions_computation(self):

        for n_points in [10, 100]:

            try:
                species_mass_fractions = stefan_tube_mass_fractions_numerical(n_points=n_points)
                (n_species, n_observations) = np.shape(species_mass_fractions)
                self.assertTrue(n_species==3)
                self.assertTrue(n_observations==n_points)

                acetone_tube_outlet = round(species_mass_fractions[0,-1], 5)
                methanol_tube_outlet = round(species_mass_fractions[1,-1], 5)
                air_tube_outlet = round(species_mass_fractions[2,-1],5)
                self.assertTrue(acetone_tube_outlet==0)
                self.assertTrue(methanol_tube_outlet==0)
                self.assertTrue(air_tube_outlet==1)

                acetone_liquid_surface = round(species_mass_fractions[0,0], 5)
                methanol_liquid_surface = round(species_mass_fractions[1,0], 5)
                self.assertTrue(acetone_liquid_surface==0.46463)
                self.assertTrue(methanol_liquid_surface==0.42424)

                idx = check.sum_of_species_fractions(species_mass_fractions, tolerance=0.000001, verbose=False)
                (idx_below_zero, idx_above_one) = check.range_of_species_fractions(species_mass_fractions, tolerance=0.000001, verbose=False)
                self.assertTrue(len(idx)==0)
                self.assertTrue(len(idx_below_zero)==0)
                self.assertTrue(len(idx_above_one)==0)
            except Exception:
                self.assertTrue(False)
