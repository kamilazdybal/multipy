import unittest
import numpy as np
import multipy

################################################################################
################################################################################
####
####    Class: Flux
####
################################################################################
################################################################################

class Flux(unittest.TestCase):

    def test_Flux__diffusive_molar_mass__allowed_calls(self):

        try:
            species_velocities = np.random.rand(2,100)
            X1 = np.random.rand(1,100)
            X2 = np.ones_like(X1) - X1
            X = np.vstack((X1, X2))
            c = np.random.rand(2,100)
            flux = multipy.Flux(species_velocities)
            df = flux.diffusive_molar_mass(X,c)
            get_df = flux.get_diffusive_molar_mass
            self.assertTrue(np.array_equal(df, get_df))
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Flux__diffusive_molar_mass__not_allowed_calls(self):

        species_velocities = np.random.rand(2,100)
        X1 = np.random.rand(1,100)
        X2 = np.ones_like(X1) - X1
        X = np.vstack((X1, X2))
        c = np.random.rand(3,100)
        flux = multipy.Flux(species_velocities)

        with self.assertRaises(ValueError):
            df = flux.diffusive_molar_mass(X,c)

        c = np.random.rand(10,2)

        with self.assertRaises(ValueError):
            df = flux.diffusive_molar_mass(X,c)

        species_velocities = np.random.rand(3,100)
        X1 = np.random.rand(1,100)
        X2 = np.ones_like(X1) - X1
        X = np.vstack((X1, X2))
        c = np.random.rand(2,100)
        flux = multipy.Flux(species_velocities)

        with self.assertRaises(ValueError):
            df = flux.diffusive_molar_mass(X,c)

################################################################################
################################################################################

    def test_Flux__diffusive_molar_mass__computation(self):

        pass

################################################################################
################################################################################
