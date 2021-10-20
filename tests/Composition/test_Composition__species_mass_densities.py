import unittest
import numpy as np
import multipy

################################################################################
################################################################################
####
####    Class: Composition
####
################################################################################
################################################################################

class Composition(unittest.TestCase):

    def test_Composition__species_mass_densities__allowed_calls(self):

        species_molar_densities = np.random.rand(10,100)
        species_molar_masses = np.random.rand(10,1)

        try:
            comp = multipy.Composition()
            rho = comp.species_mass_densities(species_molar_densities, species_molar_masses)
            (n_species, n_observations) = np.shape(rho)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 10)
        except Exception:
            self.assertTrue(False)

        species_molar_densities = np.random.rand(1,100)
        species_molar_masses = np.random.rand(1,1)

        try:
            comp = multipy.Composition()
            rho = comp.species_mass_densities(species_molar_densities, species_molar_masses)
            (n_species, n_observations) = np.shape(rho)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 1)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Composition__species_mass_densities__not_allowed_calls(self):

        species_molar_densities = np.random.rand(100,10)
        species_molar_masses = np.random.rand(10,1)
        species_molar_masses_5 = np.random.rand(5,1)
        species_molar_masses_1 = np.random.rand(1,1)
        species_molar_masses_5_0 = np.random.rand(5,1).ravel()
        species_molar_masses_1_0 = np.random.rand(1,1).ravel()
        comp = multipy.Composition()

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities(species_molar_densities, species_molar_masses_5)

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities(species_molar_densities, species_molar_masses_1)

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities(species_molar_densities, species_molar_masses_5_0)

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities(species_molar_densities, species_molar_masses_1_0)

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities(1, species_molar_masses)

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities(species_molar_densities, 1)

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities([1,2,3], species_molar_masses)

        with self.assertRaises(ValueError):
            rho = comp.species_mass_densities(species_molar_densities, [1,2,3])

################################################################################
################################################################################

    def test_Composition__species_mass_densities__computation(self):

        pass

################################################################################
################################################################################
