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

    def test_Check__range_of_species_fractions__allowed_calls(self):

        species_fractions = np.array([[0,0.1],[1,0.9]])

        try:
            check = multipy.Check()
            (idx_below_zero, idx_above_one) = check.range_of_species_fractions(species_fractions, tolerance=1e-12, verbose=False)
            self.assertTrue(np.array_equal(idx_below_zero, np.array([])))
            self.assertTrue(np.array_equal(idx_above_one, np.array([])))
        except:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Check__range_of_species_fractions__not_allowed_calls(self):

        pass

################################################################################
################################################################################
