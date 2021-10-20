import unittest
import numpy as np
import multipy

################################################################################
################################################################################
####
####    Class: Flux
####
################################################################################
################################################################################

class Flux(unittest.TestCase):

    def test__Flux__allowed_calls(self):

        species_velocities = np.random.rand(10,100)

        try:
            flux = multipy.Flux(species_velocities)
            flux.get_species_velocities
            flux.get_diffusive_molar_molar
            flux.get_diffusive_molar_mass
            flux.get_diffusive_mass_molar
            flux.get_diffusive_mass_mass
        except Exception:
            self.assertTrue(False)

        species_velocities = np.random.rand(2,100)

        try:
            flux = multipy.Flux(species_velocities)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test__Flux__not_allowed_calls(self):

        species_velocities = np.random.rand(100)

        with self.assertRaises(ValueError):
            flux = multipy.Flux(species_velocities)

        with self.assertRaises(ValueError):
            flux = multipy.Flux(1)

        with self.assertRaises(ValueError):
            flux = multipy.Flux([1,2,3])

        species_velocities = np.random.rand(1,100)

        with self.assertRaises(ValueError):
            flux = multipy.Flux(species_velocities)

################################################################################
################################################################################

    def test__Flux__allowed_set(self):

        species_velocities = np.random.rand(10,100)
        dflux = np.random.rand(10,100)

        try:
            flux = multipy.Flux(species_velocities)
            flux.set_species_velocities = species_velocities
            flux.set_diffusive_molar_molar = dflux
            flux.set_diffusive_molar_mass = dflux
            flux.set_diffusive_mass_molar = dflux
            flux.set_diffusive_mass_mass = dflux
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test__Flux__not_allowed_set(self):

        X = np.random.rand(10,100)
        flux = multipy.Flux(X)

        with self.assertRaises(AttributeError):
            flux.get_species_velocities = X

        with self.assertRaises(AttributeError):
            flux.get_diffusive_molar_molar = X

        with self.assertRaises(AttributeError):
            flux.get_diffusive_molar_mass = X

        with self.assertRaises(AttributeError):
            flux.get_diffusive_mass_molar = X

        with self.assertRaises(AttributeError):
            flux.get_diffusive_mass_mass = X

        with self.assertRaises(ValueError):
            flux.set_species_velocities = [1,2,3]

        with self.assertRaises(ValueError):
            flux.set_diffusive_molar_molar = [1,2,3]

        with self.assertRaises(ValueError):
            flux.set_diffusive_molar_mass = [1,2,3]

        with self.assertRaises(ValueError):
            flux.set_diffusive_mass_molar = [1,2,3]

        with self.assertRaises(ValueError):
            flux.set_diffusive_mass_mass = [1,2,3]

        X = np.random.rand(10,)

        with self.assertRaises(ValueError):
            flux.set_species_velocities = X

        with self.assertRaises(ValueError):
            flux.set_diffusive_molar_molar = X

        with self.assertRaises(ValueError):
            flux.set_diffusive_molar_mass = X

        with self.assertRaises(ValueError):
            flux.set_diffusive_mass_molar = X

        with self.assertRaises(ValueError):
            flux.set_diffusive_mass_mass = X

################################################################################
################################################################################
