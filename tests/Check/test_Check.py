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

    def test__Check__allowed_calls(self):

        try:
            check = multipy.Check()
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################
