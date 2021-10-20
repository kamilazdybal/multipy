import unittest
import numpy as np
import multipy

################################################################################
################################################################################
####
####    Class: Diffusion
####
################################################################################
################################################################################

class Diffusion(unittest.TestCase):

    def test__Diffusion__allowed_calls(self):

        try:
            check = multipy.Diffusion()
            self.assertTrue(check.get_n_species==0)
        except Exception:
            self.assertTrue(False)

        D1 = np.array([[1,2],
                       [2,3]])

        D2 = np.array([[1,20,10],
                       [20,3,5],
                       [10,5,4]])

        try:
            check = multipy.Diffusion(D1)
            self.assertTrue(check.get_n_species==2)
            check.set_binary_diffusion_coefficients = D2
            self.assertTrue(check.get_n_species==3)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test__Diffusion__not_allowed_calls(self):

        # Binary diffusion coefficients matrix is not symmetric:
        D1 = np.array([[1,2.1],
                       [2,3]])

        D2 = np.array([[1,10,10],
                       [20,3,5],
                       [10,5,4]])

        with self.assertRaises(ValueError):
            check = multipy.Diffusion(D1)

        with self.assertRaises(ValueError):
            check = multipy.Diffusion(D2)

        # Binary diffusion coefficients matrix is not square:
        D3 = np.array([[1,10],
                       [20,3],
                       [10,5]])

        with self.assertRaises(ValueError):
            check = multipy.Diffusion(D3)

        # Binary diffusion coefficients matrix is not numpy.ndarray
        D4 = [1,2,3]

        with self.assertRaises(ValueError):
            check = multipy.Diffusion(D4)

        D5 = 10

        with self.assertRaises(ValueError):
            check = multipy.Diffusion(D5)

################################################################################
################################################################################
