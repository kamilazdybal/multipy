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

    def test_Flux__plot_diffusive_flux__allowed_calls(self):

        try:
            species_velocities = np.random.rand(2,100)
            X1 = np.random.rand(1,100)
            X2 = np.ones_like(X1) - X1
            X = np.vstack((X1, X2))
            c = np.random.rand(2,100)
            names = [str(i) for i in range(0,2)]
            flux = multipy.Flux(species_velocities)
            df = flux.diffusive_molar_molar(X,c)
            df = flux.diffusive_molar_mass(X,c)
            df = flux.diffusive_mass_molar(X,c)
            df = flux.diffusive_mass_mass(X,c)
            flux.plot_diffusive_flux(species_names=names, figsize=(3,3), filename=None)
        except Exception:
            self.assertTrue(False)

################################################################################
################################################################################

    def test_Flux__plot_diffusive_flux__not_allowed_calls(self):

        pass

################################################################################
################################################################################
