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

    def test_Composition__plot_species_volume_fractions__allowed_calls(self):

        X = np.random.rand(2,10)

        try:
            comp = multipy.Composition()
            comp.set_species_volume_fractions = X
            plot_handle = comp.plot_species_volume_fractions()
            plot_handle.close()
        except Exception:
            self.assertTrue(False)

        try:
            comp = multipy.Composition()
            plot_handle = comp.plot_species_volume_fractions()
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Composition__plot_species_volume_fractions__not_allowed_calls(self):

        comp = multipy.Composition()
        X = np.random.rand(2,10)
        comp.set_species_volume_fractions = X

        with self.assertRaises(ValueError):
            comp.plot_species_volume_fractions(species_names=1)

        with self.assertRaises(ValueError):
            comp.plot_species_volume_fractions(custom_coordinates=1)

        with self.assertRaises(ValueError):
            comp.plot_species_volume_fractions(colors=1)

        with self.assertRaises(ValueError):
            comp.plot_species_volume_fractions(figsize=1)

        with self.assertRaises(ValueError):
            comp.plot_species_volume_fractions(filename=1)

################################################################################
################################################################################
