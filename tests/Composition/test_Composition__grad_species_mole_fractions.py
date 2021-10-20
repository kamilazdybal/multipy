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

    def test_Composition__grad_species_mole_fractions__allowed_calls(self):

        species_mole_fractions = np.random.rand(5,100)

        try:
            comp = multipy.Composition()
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=0.01)
            (n_species, n_observations) = np.shape(gradients)
            self.assertTrue(n_species==5)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

        try:
            comp = multipy.Composition()
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=1)
            (n_species, n_observations) = np.shape(gradients)
            self.assertTrue(n_species==5)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

        try:
            comp = multipy.Composition()
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=0.01, edge_order=2)
            (n_species, n_observations) = np.shape(gradients)
            self.assertTrue(n_species==5)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Composition__grad_species_mole_fractions__not_allowed_calls(self):

        species_mole_fractions = np.random.rand(5,)
        comp = multipy.Composition()

        with self.assertRaises(ValueError):
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=0.01)
        with self.assertRaises(ValueError):
            gradients = comp.grad_species_mole_fractions([1,2,3], delta=0.01)

        species_mole_fractions = np.random.rand(5,100)

        with self.assertRaises(ValueError):
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=[1])
        with self.assertRaises(ValueError):
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta='nones')

        with self.assertRaises(ValueError):
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=0.01, edge_order=0)
        with self.assertRaises(ValueError):
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=0.01, edge_order=3)
        with self.assertRaises(ValueError):
            gradients = comp.grad_species_mole_fractions(species_mole_fractions, delta=0.01, edge_order=[1])

################################################################################
################################################################################

    def test_Composition__grad_species_mole_fractions__computation(self):

        x = np.linspace(1,2,100)
        delta_x = x[2] - x[1]
        quad = x**2
        cub = x**3
        quad_grad = 2*x
        cub_grad = 3*x**2
        data = np.vstack((quad[None,:], cub[None,:]))
        expected_result = np.vstack((quad_grad[None,:], cub_grad[None,:]))

        try:
            comp = multipy.Composition()
            gradient_data = comp.grad_species_mole_fractions(data, delta_x)
            for i in range(1, 99):
                self.assertTrue(abs(gradient_data[0,i] - expected_result[0,i]) < 10e-10)
            self.assertTrue(abs(gradient_data[0,0] - expected_result[0,0]) < 0.1)
            self.assertTrue(abs(gradient_data[0,-1] - expected_result[0,-1]) < 0.1)
        except Exception:
            self.assertTrue(False)

        try:
            comp = multipy.Composition()
            gradient_data = comp.grad_species_mole_fractions(data, delta_x, edge_order=2)
            for i in range(0, 100):
                self.assertTrue(abs(gradient_data[0,i] - expected_result[0,i]) < 10e-10)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################
