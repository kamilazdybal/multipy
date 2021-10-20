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

    def test_Transform__diffusive_flux_mass_mass_to_mass_molar__allowed_calls(self):

        X = np.random.rand(5,100)
        Y = np.random.rand(5,100)

        try:
            transform = multipy.Transform()
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)
            (n_species_1,n_species_2,n_observations) = np.shape(B_uo)
            self.assertTrue(n_species_1==4)
            self.assertTrue(n_species_2==4)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

        X = np.random.rand(2,100)
        Y = np.random.rand(2,100)

        try:
            transform = multipy.Transform()
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)
            (n_species_1,n_species_2,n_observations) = np.shape(B_uo)
            self.assertTrue(n_species_1==1)
            self.assertTrue(n_species_2==1)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

        X = np.random.rand(2,1)
        Y = np.random.rand(2,1)

        try:
            transform = multipy.Transform()
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)
            (n_species_1,n_species_2,n_observations) = np.shape(B_uo)
            self.assertTrue(n_species_1==1)
            self.assertTrue(n_species_2==1)
            self.assertTrue(n_observations==1)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Transform__diffusive_flux_mass_mass_to_mass_molar__not_allowed_calls(self):

        transform = multipy.Transform()

        X = np.random.rand(1,100)
        Y = np.random.rand(1,100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)

        X = np.random.rand(5,100)
        Y = np.random.rand(4,100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)

        X = np.random.rand(5,100)
        Y = np.random.rand(1,100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)

        X = np.random.rand(1,100)
        Y = np.random.rand(5,100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)

        X = np.random.rand(100)
        Y = np.random.rand(5,100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)

        X = np.random.rand(5,100)
        Y = np.random.rand(100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, Y)

        X = np.random.rand(5,100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, [1,2,3,4,5])

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(X, None)

        Y = np.random.rand(5,100)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar([1,2,3,4,5], Y)

        with self.assertRaises(ValueError):
            B_uo = transform.diffusive_flux_mass_mass_to_mass_molar(None, Y)

################################################################################
################################################################################

    def test_Transform__diffusive_flux_mass_mass_to_mass_molar__computation(self):

        pass

################################################################################
################################################################################
