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

    def test_Velocity__plot_species_velocities__allowed_calls(self):

        species_velocities = np.random.rand(2,100)

        try:
            vel = multipy.Velocity(species_velocities)
            vel.plot_species_velocities(species_names=['1', '2'], colors=['k', 'b'], figsize=(5,5))
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Velocity__plot_species_velocities__not_allowed_calls(self):

        species_velocities = np.random.rand(2,100)
        vel = multipy.Velocity(species_velocities)

        with self.assertRaises(ValueError):
            vel.plot_species_velocities(species_names=1, colors=['k', 'b'], figsize=(5,5))

        with self.assertRaises(ValueError):
            vel.plot_species_velocities(species_names=['1', '2'], colors=1, figsize=(5,5))

        with self.assertRaises(ValueError):
            vel.plot_species_velocities(species_names=['1', '2'], colors=['k', 'b'], figsize=1)

        with self.assertRaises(ValueError):
            vel.plot_species_velocities(species_names=['1', '2', '3'], colors=['k', 'b', 'g'], figsize=(5,5))

        with self.assertRaises(ValueError):
            vel.plot_species_velocities(species_names=['1', '2'], colors=['k'], figsize=(5,5))

################################################################################
################################################################################
