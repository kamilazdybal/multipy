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

    def test_Composition__species_mole_fractions__allowed_calls(self):

        species_molar_densities = np.random.rand(10,100)
        T = 300
        p = 10000

        try:
            comp = multipy.Composition()
            X = comp.species_mole_fractions(species_molar_densities, T, p)
            (n_species, n_observations) = np.shape(X)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 10)
        except Exception:
            self.assertTrue(False)

        species_molar_densities = np.random.rand(1,100)

        try:
            comp = multipy.Composition()
            X = comp.species_mole_fractions(species_molar_densities, T, p)
            (n_species, n_observations) = np.shape(X)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 1)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Composition__species_mole_fractions__not_allowed_calls(self):

        species_molar_densities = np.random.rand(1,100).ravel()
        T = 300
        p = 10000
        comp = multipy.Composition()

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(species_molar_densities, T, p)

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(1, T, p)

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions([1,2,3], T, p)

        species_molar_densities = np.random.rand(10,100)

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(species_molar_densities, None, p)

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(species_molar_densities, 'a', p)

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(species_molar_densities, [1], p)

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(species_molar_densities, T, None)

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(species_molar_densities, T, 'a')

        with self.assertRaises(ValueError):
            X = comp.species_mole_fractions(species_molar_densities, T, [1])

################################################################################
################################################################################

    def test_Composition__species_mole_fractions__computation(self):

        species_molar_densities = np.random.rand(5,100)
        T = 10.0
        p = 8.31446261815324
        expected_result = species_molar_densities * 10.0

        try:
            comp = multipy.Composition()
            result = comp.species_mole_fractions(species_molar_densities, T, p)
            difference = abs(result - expected_result)
            residuals = difference < 1e-11
            self.assertTrue(residuals.all())
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################
