import unittest
import numpy as np
import multipy

################################################################################
################################################################################
####
####    Class: Check
####
################################################################################
################################################################################

class Check(unittest.TestCase):

    def test_Check__sum_of_species_gradients__allowed_calls(self):

        try:
            check = multipy.Check()
            X1 = np.random.rand(1,100)
            X2 = np.zeros_like(X1) - X1
            X = np.vstack((X1, X2))
            idx = np.array([])
            idx_mole = check.sum_of_species_gradients(X)
            self.assertTrue(np.array_equal(idx_mole, idx))
        except Exception:
            self.assertTrue(False)

        try:
            check = multipy.Check()
            X1 = np.array([[1, 0.4, 0.2, 0.5, 0.1]])
            X2 = np.zeros_like(X1) - np.array([[0.98, 0.35, 0.21, 0.501, 0.19]])
            X = np.vstack((X1, X2))
            idx = np.array([])
            idx_mole = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_mole, idx))
            idx_mass = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_mass, idx))
            idx_volume = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_volume, idx))
        except Exception:
            self.assertTrue(False)

        try:
            check = multipy.Check()
            X1 = np.array([[1, 0.4, 0.2, 0.5, 0.1]])
            X2 = np.zeros_like(X1) - np.array([[0.1, 0.35, 0.21, 0.501, 0.19]])
            X = np.vstack((X1, X2))
            idx = np.array([0])
            idx_mole = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_mole, idx))
            idx_mass = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_mass, idx))
            idx_volume = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_volume, idx))
        except Exception:
            self.assertTrue(False)

        try:
            check = multipy.Check()
            X1 = np.array([[1, 0.4, 0.2, 0.5, 0.1]])
            X2 = np.zeros_like(X1) - np.array([[0.1, 0.35, 0.21, 0.1, 0.19]])
            X = np.vstack((X1, X2))
            idx = np.array([0,3])
            idx_mole = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_mole, idx))
            idx_mass = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_mass, idx))
            idx_volume = check.sum_of_species_gradients(X, tolerance=0.1)
            self.assertTrue(np.array_equal(idx_volume, idx))
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Check__sum_of_species_gradients__not_allowed_calls(self):

        check = multipy.Check()
        X = np.ones((1,100))

        with self.assertRaises(ValueError):
            idx_mole = check.sum_of_species_gradients(X)

        with self.assertRaises(ValueError):
            idx_mole = check.sum_of_species_gradients([1,2,3], verbose=0)

        X1 = np.random.rand(1,100)
        X2 = np.zeros_like(X1) - X1
        X = np.vstack((X1, X2))

        with self.assertRaises(ValueError):
            idx_mole = check.sum_of_species_gradients(X, tolerance=0)

        with self.assertRaises(ValueError):
            idx_mole = check.sum_of_species_gradients(X, tolerance=1)

        with self.assertRaises(ValueError):
            idx_mole = check.sum_of_species_gradients(X, tolerance=2.5)

        with self.assertRaises(ValueError):
            idx_mole = check.sum_of_species_gradients(X, verbose=0)

################################################################################
################################################################################
