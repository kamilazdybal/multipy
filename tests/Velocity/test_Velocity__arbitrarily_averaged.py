import unittest
import numpy as np
import multipy

################################################################################
################################################################################
####
####    Class: Velocity
####
################################################################################
################################################################################

class Velocity(unittest.TestCase):

    def test_Velocity__arbitrarily_averaged__allowed_calls(self):

        species_velocities = np.random.rand(2,100)

        X1 = np.random.rand(1,100)
        X2 = np.ones_like(X1) - X1
        arbitrary_weighting_factors = np.vstack((X1, X2))

        try:
            vel = multipy.Velocity(species_velocities)
            u = vel.arbitrarily_averaged(arbitrary_weighting_factors)
            (n_dim,n_observations) = np.shape(u)
            self.assertTrue(n_dim==1)
            self.assertTrue(n_observations==100)
            u_ret = vel.get_arbitrarily_averaged
            self.assertTrue(np.array_equal(u, u_ret))
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Velocity__arbitrarily_averaged__not_allowed_calls(self):

        species_velocities = np.random.rand(2,100)

        X1 = np.random.rand(1,50)
        X2 = np.ones_like(X1) - X1
        arbitrary_weighting_factors = np.vstack((X1, X2))
        vel = multipy.Velocity(species_velocities)

        with self.assertRaises(ValueError):
            u = vel.arbitrarily_averaged(arbitrary_weighting_factors)

        X1 = np.random.rand(100)

        with self.assertRaises(ValueError):
            u = vel.arbitrarily_averaged(X1)

        with self.assertRaises(ValueError):
            vel = multipy.Velocity(1)

        with self.assertRaises(ValueError):
            vel = multipy.Velocity([1,2,3])

################################################################################
################################################################################

    def test_Velocity__arbitrarily_averaged__computation(self):

        pass

################################################################################
################################################################################
