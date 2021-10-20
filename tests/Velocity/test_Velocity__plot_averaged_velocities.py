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

    def test_Velocity__plot_averaged_velocities__allowed_calls(self):

        species_velocities = np.random.rand(2,100)
        molar_averaged = np.random.rand(1,100)
        mass_averaged = np.random.rand(1,100)
        volume_averaged = np.random.rand(1,100)
        arbitrarily_averaged = np.random.rand(1,100)

        try:
            vel = multipy.Velocity()
            vel.set_species_velocities = species_velocities
            vel.set_molar_averaged = molar_averaged
            vel.set_mass_averaged = mass_averaged
            vel.set_volume_averaged = volume_averaged
            vel.set_arbitrarily_averaged = arbitrarily_averaged
            vel.plot_averaged_velocities(colors=['k', 'b', 'g'], figsize=(10,5), filename=None)
        except Exception:
            self.assertTrue(False)

        try:
            vel = multipy.Velocity()
            vel.set_species_velocities = species_velocities
            vel.set_molar_averaged = molar_averaged
            vel.plot_averaged_velocities(colors=['k', 'b', 'g'], figsize=(10,5), filename=None)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Velocity__plot_averaged_velocities__not_allowed_calls(self):

        species_velocities = np.random.rand(2,100)
        molar_averaged = np.random.rand(1,100)
        mass_averaged = np.random.rand(1,100)
        volume_averaged = np.random.rand(1,100)
        arbitrarily_averaged = np.random.rand(1,100)
        vel = multipy.Velocity()
        vel.set_species_velocities = species_velocities
        vel.set_molar_averaged = molar_averaged
        vel.set_mass_averaged = mass_averaged
        vel.set_volume_averaged = volume_averaged
        vel.set_arbitrarily_averaged = arbitrarily_averaged

        with self.assertRaises(ValueError):
            vel.plot_averaged_velocities(colors=1, figsize=(10,5), filename=None)

        with self.assertRaises(ValueError):
            vel.plot_averaged_velocities(colors=['k', 'b', 'g'], figsize=1, filename=None)

        with self.assertRaises(ValueError):
            vel.plot_averaged_velocities(colors=['k', 'b', 'g'], figsize=(10,5), filename=1)

        with self.assertRaises(ValueError):
            vel.plot_averaged_velocities(colors=['k', 'b'], figsize=(10,5), filename=None)

################################################################################
################################################################################
