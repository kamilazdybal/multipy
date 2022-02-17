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
####    Class: Transform
####
################################################################################
################################################################################

class Transform:
    """
    Supports performing transformations of multicomponent quantities to other bases or reference frames.
    """

    def __init__(self):

        pass

    # --------------------------------------------------------------------------

    def species_fractions_mole_to_mass(self, species_mole_fractions, species_molar_masses):
        """
        Computes the species mass fractions, :math:`\\mathbf{Y}_i`, from the
        species mole fractions, :math:`\\mathbf{X}_i`, using the relation:

        .. math::

            Y_i = \\frac{M_i}{M} X_i

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions, :math:`X_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_molar_masses:
            scalar ``numpy.ndarray`` specifying the species molar masses, :math:`\\mathbf{M}_i`, in :math:`[kg/mole]`.
            It should be of size ``(n_species,1)`` where ``n_species`` is at least 2.

        :return:
            - **species_mass_fractions** - scalar ``numpy.ndarray`` specifying the species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_molar_masses, np.ndarray):
            raise ValueError("Parameter `species_molar_masses` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_dim) = np.shape(species_molar_masses)
        except:
            raise ValueError("Parameter `species_molar_masses` has to be a matrix.")

        if n_dim != 1:
            raise ValueError("Parameter `species_molar_masses` has to be of size ``(n_species,1)``.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_molar_masses` have different number of species, ``n_species``.")

        if np.any(species_molar_masses==0):
            raise ValueError("Parameter `species_molar_masses` has entries equal to zero.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mole_fractions` and `species_molar_masses` should contain all species. Only one species found.")

        composition = multipy.Composition()

        mixture_molar_mass = composition.mixture_molar_mass(species_mole_fractions, 'molar', species_molar_masses)

        species_mass_fractions = np.multiply(np.divide(species_molar_masses, mixture_molar_mass), species_mole_fractions)

        return species_mass_fractions

    # --------------------------------------------------------------------------

    def species_fractions_mass_to_mole(self, species_mass_fractions, species_molar_masses):
        """
        Computes the species mole fractions, :math:`\\mathbf{X}_i`, from the
        species mass fractions, :math:`\\mathbf{Y}_i`, using the relation:

        .. math::

            X_i = \\frac{M}{M_i} Y_i

        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mass fractions, :math:`Y_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_molar_masses:
            scalar ``numpy.ndarray`` specifying the species molar masses, :math:`\\mathbf{M}_i`, in :math:`[kg/mole]`.
            It should be of size ``(n_species,1)`` where ``n_species`` is at least 2.

        :return:
            - **species_mole_fractions** - scalar ``numpy.ndarray`` specifying the species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if not isinstance(species_molar_masses, np.ndarray):
            raise ValueError("Parameter `species_molar_masses` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_dim) = np.shape(species_molar_masses)
        except:
            raise ValueError("Parameter `species_molar_masses` has to be a matrix.")

        if n_dim != 1:
            raise ValueError("Parameter `species_molar_masses` has to be of size ``(n_species,1)``.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_molar_masses` have different number of species, ``n_species``.")

        if np.any(species_molar_masses==0):
            raise ValueError("Parameter `species_molar_masses` has entries equal to zero.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mass_fractions` and `species_molar_masses` should contain all species. Only one species found.")

        composition = multipy.Composition()

        mixture_molar_mass = composition.mixture_molar_mass(species_mass_fractions, 'mass', species_molar_masses)

        species_mole_fractions = np.multiply(np.divide(mixture_molar_mass, species_molar_masses), species_mass_fractions)

        return species_mole_fractions

    # --------------------------------------------------------------------------

    def species_gradients_mole_to_mass(self, species_mass_fractions, species_molar_masses):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`\\mathbf{J}^{XY}`,
        that allows to transform from the species mole fraction gradients, :math:`\\nabla \\mathbf{X}_i`,
        to the species mass fraction gradients, :math:`\\nabla \\mathbf{Y}_i`, according to:

        .. math::

            \\nabla \\mathbf{Y}_i = \\mathbf{J}^{XY} \\nabla \\mathbf{X}_i

        where:

        .. math::

            J_{i,j}^{XY} = \\frac{M_i}{M} \\Bigg( \\delta_{i,j} + \\frac{Y_i}{M_i} (M_n - M_j) \\Bigg)

        .. note::

            :math:`\\mathbf{J}^{XY} = (\\mathbf{J}^{YX})^{-1}`.

        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_molar_masses:
            scalar ``numpy.ndarray`` specifying **all** species molar masses, :math:`\\mathbf{M}_i`, in :math:`[kg/mole]`.
            It should be of size ``(n_species,1)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`\\mathbf{J}^{XY}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if not isinstance(species_molar_masses, np.ndarray):
            raise ValueError("Parameter `species_molar_masses` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_dim) = np.shape(species_molar_masses)
        except:
            raise ValueError("Parameter `species_molar_masses` has to be a matrix.")

        if n_dim != 1:
            raise ValueError("Parameter `species_molar_masses` has to be of size ``(n_species,1)``.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_molar_masses` have different number of species, ``n_species``.")

        if np.any(species_molar_masses==0):
            raise ValueError("Parameter `species_molar_masses` has entries equal to zero.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mass_fractions` and `species_molar_masses` should contain all species. Only one species found.")

        composition = multipy.Composition()
        mixture_molar_mass = composition.mixture_molar_mass(species_mass_fractions, 'mass', species_molar_masses)

        (n_species, n_observations) = np.shape(species_mass_fractions)

        transformation_matrix = np.zeros((n_species-1, n_species-1, n_observations))

        for k in range(0,n_observations):
            for i in range(0,n_species-1):
                for j in range(0,n_species-1):

                    if i == j:
                        kronecker_delta = 1
                    else:
                        kronecker_delta = 0

                    transformation_matrix[i,j,k] = species_molar_masses[i,0] / mixture_molar_mass[0,k] * (kronecker_delta + species_mass_fractions[i,k] / species_molar_masses[i,0] * (species_molar_masses[-1,0] - species_molar_masses[j,0]))

        return transformation_matrix

    # --------------------------------------------------------------------------

    def species_gradients_mass_to_mole(self, species_mass_fractions, species_molar_masses):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`\\mathbf{J}^{YX}`,
        that allows to transform from the species mass fraction gradients, :math:`\\nabla \\mathbf{Y}_i`,
        to the species mole fraction gradients, :math:`\\nabla \\mathbf{X}_i`, according to:

        .. math::

            \\nabla \\mathbf{X}_i = \\mathbf{J}^{YX} \\nabla \\mathbf{Y}_i

        where:

        .. math::

            J_{i,j}^{YX} = \\frac{M}{M_i} \\Bigg( \\delta_{i,j} + M Y_i \\Big( \\frac{1}{M_n} - \\frac{1}{M_j} \\Big) \\Bigg)

        .. note::

            :math:`\\mathbf{J}^{YX} = (\\mathbf{J}^{XY})^{-1}`.

        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_molar_masses:
            scalar ``numpy.ndarray`` specifying **all** species molar masses, :math:`\\mathbf{M}_i`, in :math:`[kg/mole]`.
            It should be of size ``(n_species,1)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`\\mathbf{J}^{YX}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if not isinstance(species_molar_masses, np.ndarray):
            raise ValueError("Parameter `species_molar_masses` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_dim) = np.shape(species_molar_masses)
        except:
            raise ValueError("Parameter `species_molar_masses` has to be a matrix.")

        if n_dim != 1:
            raise ValueError("Parameter `species_molar_masses` has to be of size ``(n_species,1)``.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_molar_masses` have different number of species, ``n_species``.")

        if np.any(species_molar_masses==0):
            raise ValueError("Parameter `species_molar_masses` has entries equal to zero.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mole_fractions` and `species_molar_masses` should contain all species. Only one species found.")

        composition = multipy.Composition()
        mixture_molar_mass = composition.mixture_molar_mass(species_mass_fractions, 'mass', species_molar_masses)

        (n_species, n_observations) = np.shape(species_mass_fractions)

        transformation_matrix = np.zeros((n_species-1, n_species-1, n_observations))

        for k in range(0,n_observations):
            for i in range(0,n_species-1):
                for j in range(0,n_species-1):

                    if i == j:
                        kronecker_delta = 1
                    else:
                        kronecker_delta = 0

                    transformation_matrix[i,j,k] = mixture_molar_mass[0,k] / species_molar_masses[i,0] * (kronecker_delta + mixture_molar_mass[0,k] * species_mass_fractions[i,k] * (1.0 / species_molar_masses[-1,0] - 1.0 / species_molar_masses[j,0]))

        return transformation_matrix

    # --------------------------------------------------------------------------

    def diffusive_flux_molar_molar_to_molar_volume(self, T, p, species_mole_fractions, species_partial_molar_volumes):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`\\mathbf{B}^{Vu}`,
        that allows to transform from the molar diffusive flux relative to a
        molar-averaged velocity, :math:`\\mathbf{J}_i`, to the molar diffusive flux relative
        to a volume-averaged velocity, :math:`\\mathbf{J}_i^V`, according to:

        .. math::

            \\mathbf{J}_i^V = \\mathbf{B}^{Vu} \\mathbf{J}_i

        where:

        .. math::

            B_{i,j}^{Vu} = \\delta_{i,j} - X_i (\\bar{V}_j - \\bar{V}_n) / \\bar{V}

        .. note::

            :math:`\\mathbf{B}^{Vu} = (\\mathbf{B}^{uV})^{-1}`.

        :param T: (optional)
            ``int`` or ``float`` specifying the temperature, :math:`T`, in :math:`[K]`.
        :param p: (optional)
            ``int`` or ``float`` specifying the pressure, :math:`p`, in :math:`[Pa]`.
        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_partial_molar_volumes:
            scalar ``numpy.ndarray`` specifying **all** species partial molar volumes, :math:`\\bar{\\mathbf{V}}_i`, in :math:`[m^3/mole]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`\\mathbf{B}^{Vu}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if not isinstance(T, int) and not isinstance(T, float):
                raise ValueError("Parameter `T` has to be of type `int` or `float`.")

        if not isinstance(p, int) and not isinstance(p, float):
            raise ValueError("Parameter `p` has to be of type `int` or `float`.")

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_partial_molar_volumes, np.ndarray):
            raise ValueError("Parameter `species_partial_molar_volumes` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_partial_molar_volumes)
        except:
            raise ValueError("Parameter `species_partial_molar_volumes` has to be a matrix.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_partial_molar_volumes` have different number of species, `n_species`.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_partial_molar_volumes` have different number of observations, `n_observations`.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mole_fractions` and `species_partial_molar_volumes` should contain all species. Only one species found.")

        (n_species, n_observations) = np.shape(species_mole_fractions)

        transformation_matrix = np.zeros((n_species-1, n_species-1, n_observations))

        composition = multipy.Composition()

        mixture_molar_volume = composition.mixture_molar_volume(T, p)

        for k in range(0,n_observations):
            for i in range(0,n_species-1):
                for j in range(0,n_species-1):

                    if i == j:
                        kronecker_delta = 1
                    else:
                        kronecker_delta = 0

                    transformation_matrix[i,j,k] = kronecker_delta - species_mole_fractions[i,k] * (species_partial_molar_volumes[j,k] - species_partial_molar_volumes[-1,k] ) / mixture_molar_volume

        return transformation_matrix

    # --------------------------------------------------------------------------

    def diffusive_flux_molar_volume_to_molar_molar(self, species_mole_fractions, species_partial_molar_volumes):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`\\mathbf{B}^{uV}`,
        that allows to transform from the molar diffusive flux relative to a
        volume-averaged velocity, :math:`\\mathbf{J}_i^V`, to the molar diffusive flux relative
        to a molar-averaged velocity, :math:`\\mathbf{J}_i`, according to:

        .. math::

            \\mathbf{J}_i = \\mathbf{B}^{uV} \\mathbf{J}_i^V

        where:

        .. math::

            B_{i,j}^{uV} = \\delta_{i,j} - X_i (1 - \\bar{V}_j / \\bar{V}_n)

        .. note::

            :math:`\\mathbf{B}^{uV} = (\\mathbf{B}^{Vu})^{-1}`.

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_partial_molar_volumes:
            scalar ``numpy.ndarray`` specifying **all** species partial molar volumes, :math:`\\bar{\\mathbf{V}}_i`, in :math:`[m^3/mole]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`\\mathbf{B}^{uV}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_partial_molar_volumes, np.ndarray):
            raise ValueError("Parameter `species_partial_molar_volumes` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_partial_molar_volumes)
        except:
            raise ValueError("Parameter `species_partial_molar_volumes` has to be a matrix.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_partial_molar_volumes` have different number of species, `n_species`.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_partial_molar_volumes` have different number of observations, `n_observations`.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mole_fractions` and `species_partial_molar_volumes` should contain all species. Only one species found.")

        (n_species, n_observations) = np.shape(species_mole_fractions)

        transformation_matrix = np.zeros((n_species-1, n_species-1, n_observations))

        for k in range(0,n_observations):
            for i in range(0,n_species-1):
                for j in range(0,n_species-1):

                    if i == j:
                        kronecker_delta = 1
                    else:
                        kronecker_delta = 0

                    transformation_matrix[i,j,k] = kronecker_delta - species_mole_fractions[i,k] * (1 - species_partial_molar_volumes[j,k] / species_partial_molar_volumes[-1,k])

        return transformation_matrix

    # --------------------------------------------------------------------------

    def diffusive_flux_mass_mass_to_mass_molar(self, species_mole_fractions, species_mass_fractions):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`\\mathbf{B}^{uo}`,
        that allows to transform from the mass diffusive flux relative to a
        mass-averaged velocity, :math:`\mathbf{j}_i`, to the mass diffusive flux relative
        to a molar-averaged velocity, :math:`\mathbf{j}_i^u`, according to:

        .. math::

            \\mathbf{j}_i^u = \\mathbf{B}^{uo} \\mathbf{j}_i

        where:

        .. math::

            B_{i,j}^{uo} = \\delta_{i,j} - Y_i \\Big( \\frac{X_j}{Y_j} - \\frac{X_n}{Y_n} \\Big)

        .. note::

            :math:`\\mathbf{B}^{uo} = (\\mathbf{B}^{ou})^{-1}`.

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`\\mathbf{B}^{uo}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` have different number of species, `n_species`.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` have different number of observations, `n_observations`.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` should contain all species. Only one species found.")

        (n_species, n_observations) = np.shape(species_mole_fractions)

        transformation_matrix = np.zeros((n_species-1, n_species-1, n_observations))

        for k in range(0,n_observations):
            for i in range(0,n_species-1):
                for j in range(0,n_species-1):

                    if i == j:
                        kronecker_delta = 1
                    else:
                        kronecker_delta = 0

                    transformation_matrix[i,j,k] = kronecker_delta - species_mass_fractions[i,k] * (species_mole_fractions[j,k] / species_mass_fractions[j,k] - species_mole_fractions[-1,k] / species_mass_fractions[-1,k])

        return transformation_matrix

    # --------------------------------------------------------------------------

    def diffusive_flux_mass_molar_to_mass_mass(self, species_mole_fractions, species_mass_fractions):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`\\mathbf{B}^{ou}`,
        that allows to transform from the mass diffusive flux relative to a
        molar-averaged velocity, :math:`\mathbf{j}_i^u`, to the mass diffusive flux relative
        to a mass-averaged velocity, :math:`\mathbf{j}_i`, according to:

        .. math::

            \\mathbf{j}_i = \\mathbf{B}^{ou} \\mathbf{j}_i^u

        where:

        .. math::

            B_{i,j}^{ou} = \\delta_{i,j} - Y_i \\Big( 1 - \\frac{Y_n X_j}{X_n Y_j} \\Big)

        .. note::

            :math:`\\mathbf{B}^{ou} = (\\mathbf{B}^{uo})^{-1}`.

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`\\mathbf{B}^{ou}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` have different number of species, `n_species`.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` have different number of observations, `n_observations`.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` should contain all species. Only one species found.")

        (n_species, n_observations) = np.shape(species_mole_fractions)

        transformation_matrix = np.zeros((n_species-1, n_species-1, n_observations))

        for k in range(0,n_observations):
            for i in range(0,n_species-1):
                for j in range(0,n_species-1):

                    if i == j:
                        kronecker_delta = 1
                    else:
                        kronecker_delta = 0

                    transformation_matrix[i,j,k] = kronecker_delta - species_mass_fractions[i,k] * (1 - (species_mass_fractions[-1,k] / species_mole_fractions[-1,k]) * (species_mole_fractions[j,k] / species_mass_fractions[j,k]))

        return transformation_matrix

    # --------------------------------------------------------------------------

    def fickian_diffusion_coefficients_molar_molar_to_molar_volume(self, T, p, species_mole_fractions, species_partial_molar_volumes):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`\\mathbf{B}^{Vu}`,
        that allows to transform the molar Fickian diffusion coefficients from the
        molar-averaged velocity reference frame, :math:`\mathbf{D}`, to the
        volume-averaged velocity reference frame, :math:`\mathbf{D}^V`, according to:

        .. math::

            \\mathbf{D}^V = \\mathbf{B}^{Vu} \\mathbf{D} (\\mathbf{B}^{Vu})^{-1}

        where:

        .. math::

            B_{i,j}^{Vu} = \\delta_{i,j} - X_i (\\bar{V}_j - \\bar{V}_n) / \\bar{V}

        .. note::

            :math:`\\mathbf{B}^{Vu} = (\\mathbf{B}^{uV})^{-1}`.

        :param T: (optional)
            ``int`` or ``float`` specifying the temperature, :math:`T`, in :math:`[K]`.
        :param p: (optional)
            ``int`` or ``float`` specifying the pressure, :math:`p`, in :math:`[Pa]`.
        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_partial_molar_volumes:
            scalar ``numpy.ndarray`` specifying **all** species partial molar volumes, :math:`\\bar{\\mathbf{V}}_i`, in :math:`[m^3/mole]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`\\mathbf{B}^{Vu}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        transformation_matrix = self.diffusive_flux_molar_molar_to_molar_volume(T, p, species_mole_fractions, species_partial_molar_volumes)

        return transformation_matrix

    # --------------------------------------------------------------------------

    def fickian_diffusion_coefficients_molar_molar_to_mass_mass(self, species_mole_fractions, species_mass_fractions):
        """
        Computes an invertible, :math:`n-1` dimensional transformation matrix, :math:`(\\mathbf{B}^{uo})^{-1} \\mathrm{diag}(\\mathbf{Y}_i) (\\mathrm{diag}(\\mathbf{X}_i))^{-1}`,
        that allows to transform from the molar Fickian diffusion coefficients in the
        molar-averaged velocity reference frame, :math:`\\mathbf{D}`, to the mass
        Fickian diffusion coefficients in the volume-averaged velocity reference
        frame, :math:`\\mathbf{D}^o`, according to:

        .. math::

            \\mathbf{D}^o = (\\mathbf{B}^{uo})^{-1} \\mathrm{diag}(\\mathbf{Y}_i) (\\mathrm{diag}(\\mathbf{X}_i))^{-1} \\mathbf{D} \\mathrm{diag}(\\mathbf{X}_i) (\\mathrm{diag}(\\mathbf{Y}_i))^{-1} \\mathbf{B}^{uo}

        where:

        .. math::

            B_{i,j}^{uo} = \\delta_{i,j} - Y_i \\Big( \\frac{X_j}{Y_j} - \\frac{X_n}{Y_n} \\Big)

        :math:`\\mathrm{diag}(\\mathbf{X}_i)` and :math:`\\mathrm{diag}(\\mathbf{Y}_i)` are diagonal matrices whose non-zero entries are the mole or mass fractions respectively of :math:`n-1` species.

        .. note::

            :math:`(\\mathbf{B}^{uo})^{-1} \\mathrm{diag}(\\mathbf{Y}_i) (\\mathrm{diag}(\\mathbf{X}_i))^{-1} = \\Big( \\mathrm{diag}(\\mathbf{X}_i) (\\mathrm{diag}(\\mathbf{Y}_i))^{-1} \\mathbf{B}^{uo} \\Big)^{-1}`.

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.

        :return:
            - **transformation_matrix** - scalar ``numpy.ndarray`` transformation matrix, :math:`(\\mathbf{B}^{uo})^{-1} \\mathrm{diag}(\\mathbf{Y}_i) (\\mathrm{diag}(\\mathbf{X}_i))^{-1}`, in :math:`[-]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` have different number of species, `n_species`.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` have different number of observations, `n_observations`.")

        if n_species_1 < 2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_fractions` should contain all species. Only one species found.")

        (n_species, n_observations) = np.shape(species_mole_fractions)

        transformation_matrix = np.zeros((n_species-1, n_species-1, n_observations))

        Buo = self.diffusive_flux_mass_mass_to_mass_molar(species_mole_fractions, species_mass_fractions)

        for k in range(0,n_observations):

            transformation_matrix[:,:,k] = np.dot(np.linalg.inv(Buo[:,:,k]), np.dot(np.diag(species_mass_fractions[0:-1,k].ravel()), np.linalg.inv(np.diag(species_mole_fractions[0:-1,k].ravel()))))

        return transformation_matrix
