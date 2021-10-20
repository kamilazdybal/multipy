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

    def test_Composition__species_molar_densities__allowed_calls(self):

        species_mole_fractions = np.random.rand(100,10)
        T = 300
        p = 10000

        try:
            comp = multipy.Composition()
            c = comp.species_molar_densities(species_mole_fractions, T, p)
            (n_observations, n_species) = np.shape(c)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 10)
        except Exception:
            self.assertTrue(False)

        species_mole_fractions = np.random.rand(100,1)

        try:
            comp = multipy.Composition()
            c = comp.species_molar_densities(species_mole_fractions, T, p)
            (n_observations, n_species) = np.shape(c)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 1)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Composition__species_molar_densities__not_allowed_calls(self):

        species_mole_fractions = np.random.rand(100,1).ravel()
        T = 300
        p = 10000
        comp = multipy.Composition()

        with self.assertRaises(ValueError):
            c = comp.species_molar_densities(species_mole_fractions, T, p)

        with self.assertRaises(ValueError):
            c = comp.species_molar_densities(1, T, p)

        with self.assertRaises(ValueError):
            c = comp.species_molar_densities([1,2,3], T, p)

################################################################################
################################################################################

    def test_Composition__species_molar_densities__computation(self):

        pass

################################################################################
################################################################################
