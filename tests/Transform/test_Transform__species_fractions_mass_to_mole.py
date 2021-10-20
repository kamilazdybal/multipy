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

    def test_Transform__species_fractions_mass_to_mole__allowed_calls(self):

        Y = np.random.rand(2,100)
        Mi = np.array([[1],[1]])

        try:
            transform = multipy.Transform()
            X = transform.species_fractions_mass_to_mole(Y, Mi)
            (n_species,n_observations) = np.shape(X)
            self.assertTrue(n_species==2)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

        Y = np.random.rand(5,100)
        Mi = np.array([[1],[1],[2],[2],[3]])

        try:
            transform = multipy.Transform()
            X = transform.species_fractions_mass_to_mole(Y, Mi)
            (n_species,n_observations) = np.shape(X)
            self.assertTrue(n_species==5)
            self.assertTrue(n_observations==100)
        except Exception:
            self.assertTrue(False)

        Y = np.random.rand(100,2)
        Mi = np.ones((100,1))

        try:
            transform = multipy.Transform()
            X = transform.species_fractions_mass_to_mole(Y, Mi)
            (n_species,n_observations) = np.shape(X)
            self.assertTrue(n_species==100)
            self.assertTrue(n_observations==2)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Transform__species_fractions_mass_to_mole__not_allowed_calls(self):

        Y = np.random.rand(1,100)
        Mi = np.array([[1],[1]])
        transform = multipy.Transform()

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, Mi)

        Y = np.random.rand(100,)

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, Mi)

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole([1,2], Mi)

        Y = np.random.rand(2,100)
        Mi = np.random.rand(1,2)

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, Mi)

        Mi = np.random.rand(2)

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, Mi)

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, [1,2])

        Y = np.random.rand(1,100)
        Mi = np.array([[1]])

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, Mi)

        Y = np.random.rand(2,100)
        Mi = np.array([[1],[0]])

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, Mi)

        Y = np.random.rand(2,100)
        Mi = np.array([[0],[0]])

        with self.assertRaises(ValueError):
            transform.species_fractions_mass_to_mole(Y, Mi)

################################################################################
################################################################################

    def test_Transform__species_fractions_mass_to_mole__computation_non_reacting_stefan_tube(self):

        species_mole_fractions = np.loadtxt('docs/tutorials/csv/non-reacting-stefan-tube_species-mole-fractions.csv', delimiter=',')
        species_mass_fractions = np.loadtxt('docs/tutorials/csv/non-reacting-stefan-tube_species-mass-fractions.csv', delimiter=',')
        species_molar_masses = np.array([[58.08/1000], [32.04/1000], [28.9628/1000]])
        transform = multipy.Transform()

        try:
            X_transformed = transform.species_fractions_mass_to_mole(species_mass_fractions, species_molar_masses)
            self.assertTrue(np.allclose(X_transformed, species_mole_fractions, rtol=1e-6, atol=1e-6))
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################
