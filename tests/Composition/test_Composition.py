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

    def test__Composition__allowed_calls(self):

        try:
            comp = multipy.Composition()
            comp.get_species_mole_fractions
            comp.get_species_mass_fractions
            comp.get_species_volume_fractions
            comp.get_species_molar_densities
            comp.get_species_mass_densities
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test__Composition__set_and_get(self):

        X = np.random.rand(10,100)
        Y = np.random.rand(10,100)
        Z = np.random.rand(10,100)
        C = np.random.rand(10,100)
        rho = np.random.rand(10,100)

        try:
            comp = multipy.Composition()
            comp.set_species_mole_fractions = X
            comp.set_species_mass_fractions = Y
            comp.set_species_volume_fractions = Z
            comp.set_species_molar_densities = C
            comp.set_species_mass_densities = rho
            x = comp.get_species_mole_fractions
            y = comp.get_species_mass_fractions
            z = comp.get_species_volume_fractions
            ci = comp.get_species_molar_densities
            rhoi = comp.get_species_mass_densities
            self.assertTrue(np.array_equal(x, X))
            self.assertTrue(np.array_equal(y, Y))
            self.assertTrue(np.array_equal(z, Z))
            self.assertTrue(np.array_equal(ci, C))
            self.assertTrue(np.array_equal(rhoi, rho))
            comp.set_species_mole_fractions = None
            comp.set_species_mass_fractions = None
            comp.set_species_volume_fractions = None
            comp.set_species_molar_densities = None
            comp.set_species_mass_densities = None
            x = comp.get_species_mole_fractions
            y = comp.get_species_mass_fractions
            z = comp.get_species_volume_fractions
            ci = comp.get_species_molar_densities
            rhoi = comp.get_species_mass_densities
            self.assertTrue(x==None)
            self.assertTrue(y==None)
            self.assertTrue(z==None)
            self.assertTrue(ci==None)
            self.assertTrue(rhoi==None)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test__Composition__not_allowed_set(self):

        comp = multipy.Composition()
        X = np.random.rand(10,100)

        with self.assertRaises(AttributeError):
            comp.get_species_mole_fractions = X
        with self.assertRaises(AttributeError):
            comp.get_species_mass_fractions = X
        with self.assertRaises(AttributeError):
            comp.get_species_volume_fractions = X
        with self.assertRaises(AttributeError):
            comp.get_species_molar_densities = X
        with self.assertRaises(AttributeError):
            comp.get_species_mass_densities = X

        comp = multipy.Composition()
        X = np.random.rand(10,)

        with self.assertRaises(ValueError):
            comp.set_species_mole_fractions = X
        with self.assertRaises(ValueError):
            comp.set_species_mass_fractions = X
        with self.assertRaises(ValueError):
            comp.set_species_volume_fractions = X
        with self.assertRaises(ValueError):
            comp.set_species_molar_densities = X
        with self.assertRaises(ValueError):
            comp.set_species_mass_densities = X

        with self.assertRaises(ValueError):
            comp.set_species_mole_fractions = 1
        with self.assertRaises(ValueError):
            comp.set_species_mass_fractions = 1
        with self.assertRaises(ValueError):
            comp.set_species_volume_fractions = 1
        with self.assertRaises(ValueError):
            comp.set_species_molar_densities = 1
        with self.assertRaises(ValueError):
            comp.set_species_mass_densities = 1

        with self.assertRaises(ValueError):
            comp.set_species_mole_fractions = [1,2,3]
        with self.assertRaises(ValueError):
            comp.set_species_mass_fractions = [1,2,3]
        with self.assertRaises(ValueError):
            comp.set_species_volume_fractions = [1,2,3]
        with self.assertRaises(ValueError):
            comp.set_species_molar_densities = [1,2,3]
        with self.assertRaises(ValueError):
            comp.set_species_mass_densities = [1,2,3]

################################################################################
################################################################################
