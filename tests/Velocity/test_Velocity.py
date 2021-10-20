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

    def test__Velocity__allowed_calls(self):

        species_velocities = np.random.rand(10,100)

        try:
            vel = multipy.Velocity(species_velocities)
            vel.get_species_velocities
            vel.get_molar_averaged
            vel.get_mass_averaged
            vel.get_volume_averaged
            vel.get_arbitrarily_averaged
        except Exception:
            self.assertTrue(False)

        species_velocities = np.random.rand(2,100)

        try:
            vel = multipy.Velocity(species_velocities)
        except Exception:
            self.assertTrue(False)

        try:
            vel = multipy.Velocity()
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test__Velocity__not_allowed_calls(self):

        species_velocities = np.random.rand(100)

        with self.assertRaises(ValueError):
            vel = multipy.Velocity(species_velocities)

        with self.assertRaises(ValueError):
            vel = multipy.Velocity(1)

        with self.assertRaises(ValueError):
            vel = multipy.Velocity([1,2,3])

        species_velocities = np.random.rand(1,100)

        with self.assertRaises(ValueError):
            vel = multipy.Velocity(species_velocities)

################################################################################
################################################################################

    def test__Velocity__allowed_sets(self):

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
            self.assertTrue(np.array_equal(vel.get_species_velocities, species_velocities))
            self.assertTrue(np.array_equal(vel.get_molar_averaged, molar_averaged))
            self.assertTrue(np.array_equal(vel.get_mass_averaged, mass_averaged))
            self.assertTrue(np.array_equal(vel.get_volume_averaged, volume_averaged))
            self.assertTrue(np.array_equal(vel.get_arbitrarily_averaged, arbitrarily_averaged))
        except Exception:
            self.assertTrue(False)

        try:
            vel = multipy.Velocity()
            vel.set_species_velocities = None
            vel.set_molar_averaged = None
            vel.set_mass_averaged = None
            vel.set_volume_averaged = None
            vel.set_arbitrarily_averaged = None
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test__Velocity__not_allowed_sets(self):

        species_velocities = np.random.rand(1,100)
        molar_averaged = np.random.rand(100,1)
        mass_averaged = np.random.rand(100,1)
        volume_averaged = np.random.rand(100,1)
        arbitrarily_averaged = np.random.rand(100,1)

        vel = multipy.Velocity()

        with self.assertRaises(ValueError):
            vel.set_species_velocities = species_velocities

        with self.assertRaises(ValueError):
            vel.set_molar_averaged = molar_averaged

        with self.assertRaises(ValueError):
            vel.set_mass_averaged = mass_averaged

        with self.assertRaises(ValueError):
            vel.set_volume_averaged = volume_averaged

        with self.assertRaises(ValueError):
            vel.set_arbitrarily_averaged = arbitrarily_averaged

        species_velocities = np.random.rand(100,)
        molar_averaged = np.random.rand(2,100)
        mass_averaged = np.random.rand(2,100)
        volume_averaged = np.random.rand(2,100)
        arbitrarily_averaged = np.random.rand(2,100)

        vel = multipy.Velocity()

        with self.assertRaises(ValueError):
            vel.set_species_velocities = species_velocities

        with self.assertRaises(ValueError):
            vel.set_molar_averaged = molar_averaged

        with self.assertRaises(ValueError):
            vel.set_mass_averaged = mass_averaged

        with self.assertRaises(ValueError):
            vel.set_volume_averaged = volume_averaged

        with self.assertRaises(ValueError):
            vel.set_arbitrarily_averaged = arbitrarily_averaged

        vel = multipy.Velocity()

        with self.assertRaises(ValueError):
            vel.set_species_velocities = 1

        with self.assertRaises(ValueError):
            vel.set_molar_averaged = 1

        with self.assertRaises(ValueError):
            vel.set_mass_averaged = 1

        with self.assertRaises(ValueError):
            vel.set_volume_averaged = 1

        with self.assertRaises(ValueError):
            vel.set_arbitrarily_averaged = 1

        with self.assertRaises(ValueError):
            vel.set_species_velocities = [1,2,3]

        with self.assertRaises(ValueError):
            vel.set_molar_averaged = [1,2,3]

        with self.assertRaises(ValueError):
            vel.set_mass_averaged = [1,2,3]

        with self.assertRaises(ValueError):
            vel.set_volume_averaged = [1,2,3]

        with self.assertRaises(ValueError):
            vel.set_arbitrarily_averaged = [1,2,3]

################################################################################
################################################################################
