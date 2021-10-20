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

    def test__Transform__allowed_calls(self):

        try:
            transform = multipy.Transform()
        except Exception:
            self.assertTrue(False)
