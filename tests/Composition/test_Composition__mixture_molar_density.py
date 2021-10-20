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

    def test_Composition__mixture_molar_density__allowed_calls(self):

        pass

################################################################################
################################################################################

    def test_Composition__mixture_molar_density__not_allowed_calls(self):

        pass

################################################################################
################################################################################

    def test_Composition__mixture_molar_density__computation(self):

        T = 10.0
        p = 8.31446261815324

        expected_result = 0.1

        try:
            comp = multipy.Composition()
            result = comp.mixture_molar_density(T, p)
            difference = result - expected_result
            self.assertTrue(abs(difference) < 1e-16)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################
