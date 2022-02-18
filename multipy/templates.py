"""multipy: Python library for multicomponent mass transfer"""

__author__ = "James C. Sutherland, Kamila Zdybal"
__copyright__ = "Copyright (c) 2022, James C. Sutherland, Kamila Zdybal"
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = ["Kamila Zdybal"]
__email__ = ["kamilazdybal@gmail.com"]
__status__ = "Production"

import numpy as np
import pandas as pd
import random
import copy
import scipy
import multipy
import warnings

gas_constant = 8.31446261815324

################################################################################
################################################################################
####
####    Class: Templates
####
################################################################################
################################################################################

class Templates:
    """
    Contains solution templates and building blocks for common multicomponent mass trasfer problems.
    """

    # --------------------------------------------------------------------------

    def __init__(self):

        pass

    # --------------------------------------------------------------------------

    def stefan_diffusion(self, alpha, beta):
        """
        Computes a matrix :math:`\\pmb{\\Phi}` and a vector :math:`\\pmb{\\phi}` such that:

        .. math::

            \\Phi_{i, i} = \\frac{\\alpha_i}{\\beta_{i, n}} + \\sum_{j \\neq i}^{n} \\frac{\\alpha_j}{\\beta_{i,j}}

            \\Phi_{i, j} = - \\alpha_i \\Big( \\frac{1}{\\beta_{i, j}} - \\frac{1}{\\beta_{i, n}} \\Big)

        and

        .. math::

            \\phi_i = - \\frac{\\alpha_i}{\\beta_{i, n}}

        where :math:`n` is the number of species
        and :math:`\\alpha_i` and :math:`\\beta_{i,j}` are free user-specified coefficients.

        This template can be used in Stefan diffusion -type problems.

        :param alpha:
            scalar ``numpy.ndarray`` specifying the cofficients :math:`\\alpha_{i}`. It should be of size ``(n_species,1)``.

        :param beta:
            scalar ``numpy.ndarray`` specifying the coefficients :math:`\\beta_{i,j}`. It should be of size ``(n_species,n_species)``.

        :return:
            - **phi** - scalar ``numpy.ndarray`` :math:`\\pmb{\\phi}`. It has size ``(n_species-1,1)``.
            - **Phi** - scalar ``numpy.ndarray`` :math:`\\pmb{\\Phi}`. It has size ``(n_species-1,n_species-1)``.
        """

        if not isinstance(alpha, np.ndarray):
            raise ValueError("Parameter `alpha` has to be of type `numpy.ndarray`.")

        if not isinstance(beta, np.ndarray):
            raise ValueError("Parameter `beta` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_species_2) = np.shape(beta)
        except:
            raise ValueError("Parameter `beta` has to be a matrix.")

        try:
            (n_species_3, n_dim) = np.shape(alpha)
        except:
            raise ValueError("Parameter `alpha` has to be a vector.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameter `beta` has to be a square matrix.")

        if n_species_1 != n_species_3:
            raise ValueError("Parameter `alpha` corresponds to different number of species than `beta`.")

        n_species = n_species_1

        if n_species < 2:
            raise ValueError("Parameter `beta` should correspond to all species. Only one species found.")

        Phi = np.zeros((n_species-1,n_species-1))
        phi = np.zeros((n_species-1,1))

        for i in range(0,n_species-1):
            for j in range(0,n_species-1):

                if i == j:

                    summed_terms = 0
                    for k in range(0, n_species):
                        if k != i:
                            summed_terms = summed_terms + alpha[k] / (beta[i,k])

                    Phi[i,i] = alpha[i] / (beta[i,n_species-1]) + summed_terms

                else:

                    Phi[i,j] = alpha[i] * (1/(beta[i,n_species-1]) - 1/(beta[i,j]))

            phi[i] = - alpha[i] / (beta[i,n_species-1])

        return phi, Phi

    # --------------------------------------------------------------------------
