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

    def test_Flux__diffusive_molar_molar__allowed_calls(self):

        try:
            species_velocities = np.random.rand(2,100)
            X1 = np.random.rand(1,100)
            X2 = np.ones_like(X1) - X1
            X = np.vstack((X1, X2))
            c = np.random.rand(2,100)
            flux = multipy.Flux(species_velocities)
            df = flux.diffusive_molar_molar(X,c)
            get_df = flux.get_diffusive_molar_molar
            self.assertTrue(np.array_equal(df, get_df))
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Flux__diffusive_molar_molar__not_allowed_calls(self):

        species_velocities = np.random.rand(2,100)
        X1 = np.random.rand(1,100)
        X2 = np.ones_like(X1) - X1
        X = np.vstack((X1, X2))
        c = np.random.rand(3,100)
        flux = multipy.Flux(species_velocities)

        with self.assertRaises(ValueError):
            df = flux.diffusive_molar_molar(X,c)

        c = np.random.rand(10,2)

        with self.assertRaises(ValueError):
            df = flux.diffusive_molar_molar(X,c)

        species_velocities = np.random.rand(3,100)
        X1 = np.random.rand(1,100)
        X2 = np.ones_like(X1) - X1
        X = np.vstack((X1, X2))
        c = np.random.rand(2,100)
        flux = multipy.Flux(species_velocities)

        with self.assertRaises(ValueError):
            df = flux.diffusive_molar_molar(X,c)

################################################################################
################################################################################

    def test_Flux__diffusive_molar_molar__computation(self):

        pass

################################################################################
################################################################################
