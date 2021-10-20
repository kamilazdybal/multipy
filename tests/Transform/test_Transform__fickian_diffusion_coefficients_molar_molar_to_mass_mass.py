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

    def test_Transform__fickian_diffusion_coefficients_molar_molar_to_mass_mass__allowed_calls(self):

        X = np.random.rand(5,100)
        Y = np.random.rand(5,100)

        try:
            transform = multipy.Transform()
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)
            (n_species_1,n_species_2,n_observations) = np.shape(T)
            self.assertTrue(n_species_1==4)
            self.assertTrue(n_species_2==4)
            self.assertTrue(n_observations==100)
        except:
            self.assertTrue(False)

        X = np.random.rand(2,100)
        Y = np.random.rand(2,100)

        try:
            transform = multipy.Transform()
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)
            (n_species_1,n_species_2,n_observations) = np.shape(T)
            self.assertTrue(n_species_1==1)
            self.assertTrue(n_species_2==1)
            self.assertTrue(n_observations==100)
        except:
            self.assertTrue(False)

        X = np.random.rand(2,1)
        Y = np.random.rand(2,1)

        try:
            transform = multipy.Transform()
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)
            (n_species_1,n_species_2,n_observations) = np.shape(T)
            self.assertTrue(n_species_1==1)
            self.assertTrue(n_species_2==1)
            self.assertTrue(n_observations==1)
        except:
            self.assertTrue(False)

        X = np.random.rand(100,5)
        Y = np.random.rand(100,5)

        try:
            transform = multipy.Transform()
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)
            (n_species_1,n_species_2,n_observations) = np.shape(T)
            self.assertTrue(n_species_1==99)
            self.assertTrue(n_species_2==99)
            self.assertTrue(n_observations==5)
        except:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Transform__fickian_diffusion_coefficients_molar_molar_to_mass_mass__not_allowed_calls(self):

        transform = multipy.Transform()

        X = np.random.rand(5,100)
        Y = np.random.rand(4,100)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)

        X = np.random.rand(1,100)
        Y = np.random.rand(1,100)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)

        X = np.random.rand(1,100)
        Y = np.random.rand(4,100)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)

        X = np.random.rand(4,100)
        Y = np.random.rand(1,100)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)

        X = np.random.rand(100)
        Y = np.random.rand(100)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, Y)

        X = np.random.rand(5,1)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, [1,2,3,4,5])

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(X, None)

        Y = np.random.rand(5,1)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass([1,2,3,4,5], Y)

        with self.assertRaises(ValueError):
            T = transform.fickian_diffusion_coefficients_molar_molar_to_mass_mass(None, Y)

################################################################################
################################################################################

    def test_Transform__fickian_diffusion_coefficients_molar_molar_to_mass_mass__computation(self):

        pass

################################################################################
################################################################################
