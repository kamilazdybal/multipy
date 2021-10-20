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

    def test_Composition__species_mass_fractions__allowed_calls(self):

        species_mass_densities = np.random.rand(10,100)
        mixture_mass_density = np.random.rand(1,100)

        try:
            comp = multipy.Composition()
            Y = comp.species_mass_fractions(species_mass_densities, mixture_mass_density)
            (n_species, n_observations) = np.shape(Y)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 10)
        except Exception:
            self.assertTrue(False)

        species_mass_densities = np.random.rand(1,100)
        mixture_mass_density = np.random.rand(1,100)

        try:
            comp = multipy.Composition()
            Y = comp.species_mass_fractions(species_mass_densities, mixture_mass_density)
            (n_species, n_observations) = np.shape(Y)
            self.assertTrue(n_observations == 100)
            self.assertTrue(n_species == 1)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Composition__species_mass_fractions__not_allowed_calls(self):

        species_mass_densities = np.random.rand(10,100)
        species_mass_densities_0 = np.random.rand(100,)
        mixture_mass_density = np.random.rand(1,100)
        mixture_mass_density_0 = np.random.rand(100,)
        mixture_mass_density_2 = np.random.rand(2,100)
        comp = multipy.Composition()

        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(species_mass_densities_0, mixture_mass_density)
        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(species_mass_densities[:,0:10], mixture_mass_density)
        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(1, mixture_mass_density)
        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions([1,2,3], mixture_mass_density)

        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(species_mass_densities, mixture_mass_density_0)
        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(species_mass_densities, mixture_mass_density_2)
        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(species_mass_densities, mixture_mass_density[:,0:10])
        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(species_mass_densities, 1)
        with self.assertRaises(ValueError):
            Y = comp.species_mass_fractions(species_mass_densities, [1,2,3])

################################################################################
################################################################################

    def test_Composition__species_mass_fractions__computation(self):

        species_mass_densities = np.array([[2,20,200],[1,10,100]])
        mixture_mass_density = np.array([[1,10,100]])
        expected_result = np.array([[2,2,2],[1,1,1]])

        try:
            comp = multipy.Composition()
            result = comp.species_mass_fractions(species_mass_densities, mixture_mass_density)
            difference = abs(result - expected_result)
            residuals = difference < 1e-11
            self.assertTrue(residuals.all())
        except Exception:
            self.assertTrue(False)

        species_mass_densities = np.array([[2,20],[1,10]])
        mixture_mass_density = np.array([[1000,10]])
        expected_result = np.array([[0.002, 2],[0.001,1]])

        try:
            comp = multipy.Composition()
            result = comp.species_mass_fractions(species_mass_densities, mixture_mass_density)
            difference = abs(result - expected_result)
            residuals = difference < 1e-11
            self.assertTrue(residuals.all())
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################
