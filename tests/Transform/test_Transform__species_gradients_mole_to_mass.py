import unittest
import numpy as np
import multipy

################################################################################
################################################################################
####
####    Class: Transform
####
################################################################################
################################################################################

class Transform(unittest.TestCase):

    def test_Transform__species_gradients_mole_to_mass__allowed_calls(self):

        grad_X = np.random.rand(2,100)
        Mi = np.array([[1],[1]])

        try:
            transform = multipy.Transform()
            J_XY = transform.species_gradients_mole_to_mass(grad_X, Mi)
            (n_species_1,n_species_2,n_observations) = np.shape(J_XY)
            self.assertTrue(n_species_1==1)
            self.assertTrue(n_species_2==1)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

        grad_X = np.random.rand(5,100)
        Mi = np.array([[1],[1],[2],[2],[3]])

        try:
            transform = multipy.Transform()
            J_XY = transform.species_gradients_mole_to_mass(grad_X, Mi)
            (n_species_1,n_species_2,n_observations) = np.shape(J_XY)
            self.assertTrue(n_species_1==4)
            self.assertTrue(n_species_2==4)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Transform__species_gradients_mole_to_mass__not_allowed_calls(self):

        transform = multipy.Transform()

        grad_X = np.random.rand(1,100)
        Mi = np.array([[1]])

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass(grad_X, Mi)

        grad_X = np.random.rand(5,100)
        Mi = np.array([[1],[2],[3]])

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass(grad_X, Mi)

        grad_X = np.random.rand(5,100)
        Mi = np.array([1,2,3,4,5])

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass(grad_X, Mi)

        grad_X = np.random.rand(5)
        Mi = np.array([[1],[1],[2],[2],[3]])

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass(grad_X, Mi)

        grad_X = np.random.rand(5,100)

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass(grad_X, [1,2,3,4,5])

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass(grad_X, None)

        Mi = np.array([[1],[1],[2],[2],[3]])

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass([1,2,3,4,5], Mi)

        with self.assertRaises(ValueError):
            J_XY = transform.species_gradients_mole_to_mass(None, Mi)

################################################################################
################################################################################

    def test_Transform__species_gradients_mole_to_mass__computation_non_reacting_stefan_tube(self):

        species_mole_fractions = np.loadtxt('docs/tutorials/csv/non-reacting-stefan-tube_species-mole-fractions.csv', delimiter=',')
        species_mass_fractions = np.loadtxt('docs/tutorials/csv/non-reacting-stefan-tube_species-mass-fractions.csv', delimiter=',')
        species_molar_masses = np.array([[58.08/1000], [32.04/1000], [28.9628/1000]])
        transform = multipy.Transform()
        composition = multipy.Composition()
        delta = 0.001001
        grad_species_mole_fractions = composition.grad_species_mole_fractions(species_mole_fractions, delta, edge_order=2)
        grad_species_mass_fractions = composition.grad_species_mass_fractions(species_mass_fractions, delta, edge_order=2)

        (n_species, n_observations) = np.shape(species_mole_fractions)

        try:
            JXY = transform.species_gradients_mole_to_mass(species_mass_fractions, species_molar_masses)
            grad_species_mass_fractions_transformed = np.zeros_like(grad_species_mass_fractions)
            for k in range(0,n_observations):
                grad_species_mass_fractions_transformed[0:n_species-1,k] = np.dot(JXY[:,:,k], grad_species_mole_fractions[0:n_species-1,k])

            for k in range(0,n_observations):
                grad_species_mass_fractions_transformed[-1,k] = - np.sum(grad_species_mass_fractions_transformed[0:-1,k])

            self.assertTrue(np.allclose(grad_species_mass_fractions, grad_species_mass_fractions_transformed, rtol=1e-5, atol=1e-5))
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################
