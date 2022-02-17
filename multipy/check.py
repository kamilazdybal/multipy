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
####    Class: Check
####
################################################################################
################################################################################

class Check:
    """
    Supports performing basic checks of the computed quantities.
    """

    # --------------------------------------------------------------------------

    def __init__(self):

        pass

    # --------------------------------------------------------------------------

    def sum_of_species_fractions(self, species_fractions, tolerance=1e-12, verbose=False):
        """
        Checks if all species mole/mass/volume fractions sum to 1.0 for
        every observation within a specified tolerance.

        For mole fractions:

        .. math::

            \\sum_{i=1}^{n} X_i = 1.0

        For mass fractions:

        .. math::

            \\sum_{i=1}^{n} Y_i = 1.0

        For volume fractions:

        .. math::

            \\sum_{i=1}^{n} V_i = 1.0

        where :math:`n` is the number of species.

        :param species_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole/mass/volume fractions in :math:`[-]`.
            It should be of size ``(n_species, n_observations)`` where
            ``n_species`` is at least 2.
        :param tolerance: (optional)
            ``float`` specifying the tolerance. It should be larger than 0.0 and
            smaller than 1.0.
        :param verbose: (optional)
            ``bool`` for printing verbose information.

        :return:
            - **idx** - indices of observations where species mole/mass/volume fractions do not sum to 1.0 within a specified tolerance.
        """

        if not isinstance(species_fractions, np.ndarray):
            raise ValueError("Parameter `species_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_fractions)
        except:
            raise ValueError("Parameter `species_fractions` has to be a matrix.")

        if n_species < 2:
            raise ValueError("Species fractions matrix `species_mole_fractions` has to have at least two species.")

        if n_observations < n_species:
            warnings.warn("Number of observations in `species_fractions` is smaller than the number of species. Make sure that the `species_fractions` has shape `(n_observations,n_species)`.")

        if not isinstance(tolerance, float):
            raise ValueError("Parameter `tolerance` has to be of type `float`.")

        if tolerance <= 0 or tolerance >= 1:
            raise ValueError("Parameter `tolerance` has to be larger than 0 and smaller than 1.")

        if not isinstance(verbose, bool):
            raise ValueError("Parameter `verbose` has to be of type `bool`.")

        sums = np.sum(species_fractions, axis=0)
        sums_boolean = np.zeros_like(sums)

        for i, observation in enumerate(sums):

            if (observation < 1+tolerance) and (observation > 1-tolerance):
                sums_boolean[i] = True
            else:
                sums_boolean[i] = False

        if sums_boolean.all():
            if verbose: print('All mole/mass/volume fractions sum to 1.0 within a specified tolerance.')
            idx = np.array([])
        else:
            if verbose: print('Detected observations where mole/mass/volume fractions do not sum to 1.0 within a specified tolerance.')
            (idx, ) = np.where(sums_boolean==False)

        return idx

    # --------------------------------------------------------------------------

    def range_of_species_fractions(self, species_fractions, tolerance=1e-12, verbose=False):
        """
        Checks if all species mole/mass/volume fraction values are bounded between
        0 and 1.

        For mole fractions:

        .. math::

            X_i \\in \\langle 0, 1 \\rangle

        For mass fractions:

        .. math::

            Y_i \\in \\langle 0, 1 \\rangle

        For volume fractions:

        .. math::

            V_i \\in \\langle 0, 1 \\rangle

        :param species_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole/mass/volume fractions in :math:`[-]`.
            It should be of size ``(n_observations,n_species)`` where
            ``n_species`` is at least 2.
        :param verbose: (optional)
            ``bool`` for printing verbose information.

        :return:
            - **idx_below_zero** - indices of observations where species mole/mass/volume fractions are less than 0.0 within a specified tolerance.
            - **idx_above_one** - indices of observations where species mole/mass/volume fractions are larger than 1.0 within a specified tolerance.
        """

        if not isinstance(species_fractions, np.ndarray):
            raise ValueError("Parameter `species_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_fractions)
        except:
            raise ValueError("Parameter `species_fractions` has to be a matrix.")

        if n_species < 2:
            raise ValueError("Mole fractions matrix `species_fractions` has to have at least two species.")

        if n_observations < n_species:
            warnings.warn("Number of observations in `species_fractions` is smaller than the number of species. Make sure that the `species_fractions` has shape `(n_observations,n_species)`.")

        if not isinstance(verbose, bool):
            raise ValueError("Parameter `verbose` has to be of type `bool`.")

        if not np.greater_equal(species_fractions, 0-tolerance).all():
            if verbose: print('Not all mole/mass/volume fractions are larger than 0.0 within a specified tolerance.')
            (idx_below_zero_i, idx_below_zero_j) = np.where(species_fractions<(0-tolerance))
            idx_below_zero = np.hstack((idx_below_zero_i[:,None], idx_below_zero_j[:,None]))
        else:
            if verbose: print('All mole/mass/volume fractions are larger than 0.0 within a specified tolerance.')
            idx_below_zero = np.array([])

        if not np.less_equal(species_fractions, 1+tolerance).all():
            if verbose: print('Not all mole/mass/volume fractions are smaller than 1.0 within a specified tolerance.')
            (idx_above_one_i, idx_above_one_j) = np.where(species_fractions>(1+tolerance))
            idx_above_one = np.hstack((idx_above_one_i[:,None], idx_above_one_j[:,None]))
        else:
            if verbose: print('All mole/mass/volume fractions are smaller than 1.0 within a specified tolerance.')
            idx_above_one = np.array([])

        return (idx_below_zero, idx_above_one)

    # --------------------------------------------------------------------------

    def sum_of_species_gradients(self, species_gradients, tolerance=1e-12, verbose=False):
        """
        Checks if all species mole/mass/volume fraction gradients sum to 0.0 for
        every observation within a specified tolerance.

        For mole fractions:

        .. math::

            \\sum_{i=1}^{n} \\nabla X_i = 0.0

        For mass fractions:

        .. math::

            \\sum_{i=1}^{n} \\nabla Y_i = 0.0

        For volume fractions:

        .. math::

            \\sum_{i=1}^{n} \\nabla V_i = 0.0

        where :math:`n` is the number of species.

        :param species_gradients:
            scalar ``numpy.ndarray`` specifying **all** species mole/mass/volume fraction gradients in :math:`[-]`.
            It should be of size ``(n_species, n_observations)`` where
            ``n_species`` is at least 2.
        :param tolerance: (optional)
            ``float`` specifying the tolerance. It should be larger than 0.0 and
            smaller than 1.0.
        :param verbose: (optional)
            ``bool`` for printing verbose information.

        :return:
            - **idx** - indices of observations where species mole/mass/volume fraction gradients do not sum to 0.0 within a specified tolerance.
        """

        if not isinstance(species_gradients, np.ndarray):
            raise ValueError("Parameter `species_gradients` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_gradients)
        except:
            raise ValueError("Parameter `species_gradients` has to be a matrix.")

        if n_species < 2:
            raise ValueError("Species fractions matrix `species_gradients` has to have at least two species.")

        if n_observations < n_species:
            warnings.warn("Number of observations in `species_gradients` is smaller than the number of species. Make sure that the `species_fractions` has shape `(n_observations,n_species)`.")

        if not isinstance(tolerance, float):
            raise ValueError("Parameter `tolerance` has to be of type `float`.")

        if tolerance <= 0 or tolerance >= 1:
            raise ValueError("Parameter `tolerance` has to be larger than 0 and smaller than 1.")

        if not isinstance(verbose, bool):
            raise ValueError("Parameter `verbose` has to be of type `bool`.")

        sums = np.sum(species_gradients, axis=0)
        sums_boolean = np.zeros_like(sums)

        for i, observation in enumerate(sums):

            if (observation < tolerance) and (observation > -tolerance):
                sums_boolean[i] = True
            else:
                sums_boolean[i] = False

        if sums_boolean.all():
            if verbose: print('All mole/mass/volume fraction gradiens sum to 0.0 within a specified tolerance.')
            idx = np.array([])
        else:
            if verbose: print('Detected observations where mole/mass/volume fraction gradients do not sum to 0.0 within a specified tolerance.')
            (idx, ) = np.where(sums_boolean==False)

        return idx

    # --------------------------------------------------------------------------

    def sum_of_species_production_rates(self, species_production_rates, tolerance=1e-12, verbose=False):
        """
        Checks if all species production rates sum to 0.0 for
        every observation within a specified tolerance:

        For net molar production rates:

        .. math::

            \\sum_{i=1}^{n} s_i = 0.0

        For net mass production rates:

        .. math::

            \\sum_{i=1}^{n} \\omega_i = 0.0

        where :math:`n` is the number of species.

        :param species_production_rates:
            scalar ``numpy.ndarray`` specifying **all** species production rates, :math:`s_i` in :math:`mole/(m^3s)` or :math:`\\omega_i` in :math:`kg/(m^3s)`.
            It should be of size ``(n_species,n_observations)`` where
            ``n_species`` is at least 2.
        :param tolerance: (optional)
            ``float`` specifying the tolerance. It should be larger than 0.0 and
            smaller than 1.0.
        :param verbose: (optional)
            ``bool`` for printing verbose information.

        :return:
            - **idx** - indices of observations where species source terms do not sum to 0.0 within a specified tolerance.
        """

        if not isinstance(species_production_rates, np.ndarray):
            raise ValueError("Parameter `species_production_rates` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_production_rates)
        except:
            raise ValueError("Parameter `species_production_rates` has to be a matrix.")

        if n_species < 2:
            raise ValueError("Species source terms matrix `species_production_rates` has to have at least two species.")

        if n_observations < n_species:
            warnings.warn("Number of observations in `species_production_rates` is smaller than the number of species. Make sure that the `species_production_rates` has shape `(n_observations,n_species)`.")

        if not isinstance(tolerance, float):
            raise ValueError("Parameter `tolerance` has to be of type `float`.")

        if tolerance <= 0 or tolerance >= 1:
            raise ValueError("Parameter `tolerance` has to be larger than 0 and smaller than 1.")

        if not isinstance(verbose, bool):
            raise ValueError("Parameter `verbose` has to be of type `bool`.")

        sums = np.sum(species_production_rates, axis=0)
        sums_boolean = np.zeros_like(sums)

        for i, observation in enumerate(sums):

            if (observation < tolerance) and (observation > -tolerance):
                sums_boolean[i] = True
            else:
                sums_boolean[i] = False

        if sums_boolean.all():
            if verbose: print('All species production rates sum to 0.0 within a specified tolerance.')
            idx = np.array([])
        else:
            if verbose: print('Detected observations where species production rates do not sum to 0.0 within a specified tolerance.')
            (idx, ) = np.where(sums_boolean==False)

        return idx

    # --------------------------------------------------------------------------
