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
####    Class: Diffusion
####
################################################################################
################################################################################

class Diffusion:
    """
    Supports computing and storing diffusion-related quantities.

    :param binary_diffusion_coefficients:
        scalar ``numpy.ndarray`` specifying the binary diffusion coefficients, :math:`\\pmb{\\mathcal{D}}`, in :math:`[m^2/s]` for **all** species. It should be a symmetric matrix of size ``(n_species,n_species)``.

    **Getters:**

    - **get_binary_diffusion_coefficients** (is set to ``None`` at class init)
    - **get_n_species** read only - (is set to ``0`` at class init)
    - **get_fickian_diffusion_coefficients_molar_molar** (is set to ``None`` at class init)
    - **get_fickian_diffusion_coefficients_mass_mass** (is set to ``None`` at class init)
    - **get_fickian_diffusion_coefficients_molar_volume** (is set to ``None`` at class init)

    **Setters:**

    - **set_binary_diffusion_coefficients** setter for ``get_binary_diffusion_coefficients``
    - **set_fickian_diffusion_coefficients_molar_molar** setter for ``get_fickian_diffusion_coefficients_molar_molar``
    - **set_fickian_diffusion_coefficients_mass_mass** setter for ``get_fickian_diffusion_coefficients_mass_mass``
    - **set_fickian_diffusion_coefficients_molar_volume** setter for ``get_fickian_diffusion_coefficients_molar_volume``
    """

    # --------------------------------------------------------------------------

    def __init__(self, binary_diffusion_coefficients=None):

        if binary_diffusion_coefficients is not None:
            if not isinstance(binary_diffusion_coefficients, np.ndarray):
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be of type `numpy.ndarray`.")

            try:
                (n_coefficients_1, n_coefficients_2) = np.shape(binary_diffusion_coefficients)
            except:
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be a matrix.")

            if n_coefficients_1 != n_coefficients_2:
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be a square matrix.")

            if not np.allclose(binary_diffusion_coefficients, binary_diffusion_coefficients.T, rtol=1e-08, atol=1e-08):
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be a symmetric matrix.")

            self.__n_species = n_coefficients_1

        else:
            self.__n_species = 0

        self.__binary_diffusion_coefficients = binary_diffusion_coefficients
        self.__fickian_diffusion_coefficients_molar_molar = None
        self.__fickian_diffusion_coefficients_mass_mass = None
        self.__fickian_diffusion_coefficients_molar_volume = None

    @property
    def get_binary_diffusion_coefficients(self):
        return self.__binary_diffusion_coefficients

    @property
    def get_n_species(self):
        return self.__n_species

    @property
    def get_fickian_diffusion_coefficients_molar_molar(self):
        return self.__fickian_diffusion_coefficients_molar_molar

    @property
    def get_fickian_diffusion_coefficients_mass_mass(self):
        return self.__fickian_diffusion_coefficients_mass_mass

    @property
    def get_fickian_diffusion_coefficients_molar_volume(self):
        return self.__fickian_diffusion_coefficients_molar_volume

    @get_binary_diffusion_coefficients.setter
    def set_binary_diffusion_coefficients(self, new_binary_diffusion_coefficients):

        if new_binary_diffusion_coefficients is not None:
            if not isinstance(new_binary_diffusion_coefficients, np.ndarray):
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be of type `numpy.ndarray`.")

            try:
                (n_coefficients_1, n_coefficients_2) = np.shape(new_binary_diffusion_coefficients)
            except:
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be a matrix.")

            if n_coefficients_1 != n_coefficients_2:
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be a square matrix.")

            if not np.allclose(new_binary_diffusion_coefficients, new_binary_diffusion_coefficients.T, rtol=1e-08, atol=1e-08):
                raise ValueError("Parameter `binary_diffusion_coefficients` has to be a symmetric matrix.")

            self.__n_species = n_coefficients_1

        self.__binary_diffusion_coefficients = new_binary_diffusion_coefficients

    @get_fickian_diffusion_coefficients_molar_molar.setter
    def set_fickian_diffusion_coefficients_molar_molar(self, new_fickian_diffusion_coefficients):

        if new_fickian_diffusion_coefficients is not None:
            if not isinstance(new_fickian_diffusion_coefficients, np.ndarray):
                raise ValueError("Parameter `fickian_diffusion_coefficients_molar_molar` has to be of type `numpy.ndarray`.")

            try:
                (n_coefficients_1, n_coefficients_2, n_observations) = np.shape(new_fickian_diffusion_coefficients)
            except:
                raise ValueError("Parameter `fickian_diffusion_coefficients_molar_molar` has to be a 3D matrix.")

        if n_coefficients_1 != n_coefficients_2:
            raise ValueError("Parameter `fickian_diffusion_coefficients_molar_molar` has to be a square matrix.")

        self.__fickian_diffusion_coefficients_molar_molar = new_fickian_diffusion_coefficients

    @get_fickian_diffusion_coefficients_mass_mass.setter
    def set_fickian_diffusion_coefficients_mass_mass(self, new_fickian_diffusion_coefficients):

        if new_fickian_diffusion_coefficients is not None:
            if not isinstance(new_fickian_diffusion_coefficients, np.ndarray):
                raise ValueError("Parameter `fickian_diffusion_coefficients_mass_mass` has to be of type `numpy.ndarray`.")

            try:
                (n_coefficients_1, n_coefficients_2, n_observations) = np.shape(new_fickian_diffusion_coefficients)
            except:
                raise ValueError("Parameter `fickian_diffusion_coefficients_mass_mass` has to be a 3D matrix.")

        if n_coefficients_1 != n_coefficients_2:
            raise ValueError("Parameter `fickian_diffusion_coefficients_mass_mass` has to be a square matrix.")

        self.__fickian_diffusion_coefficients_mass_mass = new_fickian_diffusion_coefficients

    @get_fickian_diffusion_coefficients_molar_volume.setter
    def set_fickian_diffusion_coefficients_molar_volume(self, new_fickian_diffusion_coefficients):

        if new_fickian_diffusion_coefficients is not None:
            if not isinstance(new_fickian_diffusion_coefficients, np.ndarray):
                raise ValueError("Parameter `fickian_diffusion_coefficients_molar_volume` has to be of type `numpy.ndarray`.")

            try:
                (n_coefficients_1, n_coefficients_2, n_observations) = np.shape(new_fickian_diffusion_coefficients)
            except:
                raise ValueError("Parameter `fickian_diffusion_coefficients_molar_volume` has to be a 3D matrix.")

        if n_coefficients_1 != n_coefficients_2:
            raise ValueError("Parameter `fickian_diffusion_coefficients_molar_volumer` has to be a square matrix.")

        self.__fickian_diffusion_coefficients_molar_volume = new_fickian_diffusion_coefficients

    # --------------------------------------------------------------------------

    def fickian_diffusion_coefficients_molar_molar(self, species_mole_fractions):
        """
        Computes the molar Fickian diffusion coefficients expressed in a molar-averaged velocity reference frame, :math:`\mathbf{D}`, from:

        .. math::

            \\mathbf{D} = \\mathbf{B}^{-1}

        where the elements of :math:`\\mathbf{B}` are given by:

        .. math::

            B_{i, i} = \\frac{X_i}{\\mathcal{D}_{i, n}} + \\sum_{j \\neq i}^{n} \\frac{X_j}{\\mathcal{D}_{i,j}}

            B_{i, j} = - X_i \\Big( \\frac{1}{\\mathcal{D}_{i, j}} - \\frac{1}{\\mathcal{D}_{i, n}} \\Big)

        and :math:`n` is the number of species.

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying the species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.

        :return:
            - **fickian_diffusion_coefficients** - scalar ``numpy.ndarray`` of molar Fickian diffusion coefficients expressed in a molar-averaged velocity reference frame, :math:`\\mathbf{D}`, in :math:`[m^2/s]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        if self.get_binary_diffusion_coefficients is None:
            raise ValueError("Parameter `binary_diffusion_coefficients` has to be set in the `Diffusion` class object.")

        (n_species, n_observations) = np.shape(species_mole_fractions)

        if n_species != self.get_n_species or n_species != self.get_n_species:
            raise ValueError("Parameters `species_mole_fractions` and `binary_diffusion_coefficients` have different number of species `n_species`.")

        fickian_diffusion_coefficients = np.zeros((self.get_n_species-1, self.get_n_species-1, n_observations))

        for observation in range(0,n_observations):

            current_composition = species_mole_fractions[:,observation]

            B = np.zeros((self.get_n_species-1, self.get_n_species-1))

            for i in range(0,self.get_n_species-1):
                for j in range(0,self.get_n_species-1):
                    if i == j:
                        summed_terms = 0
                        for k in range(0, self.get_n_species):
                            if k != i:
                                summed_terms = summed_terms + current_composition[k] / self.get_binary_diffusion_coefficients[i, k]
                        B[i,j] = current_composition[i] / self.get_binary_diffusion_coefficients[i,self.get_n_species-1] + summed_terms
                    else:
                        B[i,j] = - current_composition[i] * (1/self.get_binary_diffusion_coefficients[i,j] - 1/self.get_binary_diffusion_coefficients[i,self.get_n_species-1])

            current_fickian_diffusion_coefficients = np.linalg.inv(B)

            fickian_diffusion_coefficients[:,:,observation] = current_fickian_diffusion_coefficients

        self.__fickian_diffusion_coefficients_molar_molar = fickian_diffusion_coefficients

        return fickian_diffusion_coefficients

    # --------------------------------------------------------------------------

    def fickian_diffusion_coefficients_mass_mass(self, species_mole_fractions, species_mass_fractions):
        """
        Computes the mass Fickian diffusion coefficients expressed in a mass-averaged velocity reference frame, :math:`\\mathbf{D}^o`, from:

        .. math::

            \\mathbf{D}^o =

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying the species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying the species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.

        :return:
            - **fickian_diffusion_coefficients** - scalar ``numpy.ndarray`` of mass Fickian diffusion expressed in a mass-averaged velocity reference frame, :math:`\\mathbf{D}^o`, in :math:`[m^2/s]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        pass

    # --------------------------------------------------------------------------

    def fickian_diffusion_coefficients_molar_volume(self, species_volume_fractions):
        """
        Computes the molar Fickian diffusion coefficients expressed in a volume-averaged velocity reference frame, :math:`\mathbf{D}^V`, from:

        .. math::

            \\mathbf{D}^V =

        :param species_volume_fractions:
            scalar ``numpy.ndarray`` specifying the species volume fractions, :math:`\\mathbf{V}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.

        :return:
            - **fickian_diffusion_coefficients** - scalar ``numpy.ndarray`` of molar Fickian diffusion expressed in a volume-averaged velocity reference frame, :math:`\\mathbf{D}^V`, in :math:`[m^2/s]`. It has size ``(n_species-1,n_species-1,n_observations)``.
        """

        pass

    # --------------------------------------------------------------------------

    def effective_diffusivity(self, case=None):
        """
        Computes the effective diffusivity using one of the selected assumptions:


        """

        pass

    # --------------------------------------------------------------------------

    def diffusive_flux_molar_molar(self, grad_species_mole_fractions, mixture_molar_density):
        """
        Computes the molar diffusive flux relative to a molar-averaged velocity, :math:`\\mathbf{J}_i`,
        from the generalized Fick's law:

        .. math::

            \\mathbf{J}_i = - c \\mathbf{D} \\nabla \\mathbf{X}_i

        where :math:`\\mathbf{D}` are the molar Fickian diffusion coefficients expressed in a molar-averaged velocity reference frame.

        :param grad_species_mole_fractions:
            vector ``numpy.ndarray`` specifying the species mole fractions gradients, :math:`\\nabla \\mathbf{X}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param mixture_molar_density:
            mixture molar density, :math:`c`, in :math:`[mole/m^3]`.

        :return:
            - **diffusive_flux** - vector ``numpy.ndarray`` of molar diffusive fluxes relative to a molar-averaged velocity, :math:`\\mathbf{J}_i`, in :math:`[mole/(m^2s)]`. It has size ``(n_species,n_observations)``.
        """

        if self.get_fickian_diffusion_coefficients_molar_molar is None:
            raise ValueError("Parameter `fickian_diffusion_coefficients_molar_molar` has to be set in the `Diffusion` class object.")

        (n_species, n_observations) = np.shape(grad_species_mole_fractions)

        diffusive_flux_independent = np.zeros((self.get_n_species-1, n_observations))

        for index in range(0,n_observations):

            diffusive_flux_independent[:,index] = - mixture_molar_density * np.dot(self.get_fickian_diffusion_coefficients_molar_molar[:,:,index], grad_species_mole_fractions[0:-1,index])

        diffusive_flux_last_species = - np.sum(diffusive_flux_independent, axis=0)

        diffusive_flux = np.vstack((diffusive_flux_independent, diffusive_flux_last_species[None,:]))

        return diffusive_flux

    # --------------------------------------------------------------------------

    def diffusive_flux_mass_mass(self, grad_species_mass_fractions, mixture_mass_density):
        """
        Computes the mass diffusive flux relative to a mass-averaged velocity, :math:`\\mathbf{j}_i`,
        from the generalized Fick's law:

        .. math::

            \\mathbf{j}_i = - \\rho \\mathbf{D}^o \\nabla \\mathbf{Y}_i

        where :math:`\\mathbf{D}^o` are the mass Fickian diffusion coefficients expressed in a molar-averaged velocity reference frame.

        :param grad_species_mass_fractions:
            vector ``numpy.ndarray`` specifying the species mass fractions gradients, :math:`\\nabla \\mathbf{Y}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param mixture_mass_density:
            scalar ``numpy.ndarray`` specifying the mixture mass density, :math:`\\pmb{\\rho}`, in :math:`[kg/m^3]`.
            It should be of size ``(1,n_observations)``.

        :return:
            - **diffusive_flux** - vector ``numpy.ndarray`` of mass diffusive fluxes relative to a mass-averaged velocity, :math:`\\mathbf{j}_i`, in :math:`[kg/(m^2s)]`. It has size ``(n_species,n_observations)``.
        """

        if self.get_fickian_diffusion_coefficients_mass_mass is None:
            raise ValueError("Parameter `fickian_diffusion_coefficients_mass_mass` has to be set in the `Diffusion` class object.")

        (n_species, n_observations) = np.shape(grad_species_mass_fractions)

        diffusive_flux_independent = np.zeros((self.get_n_species-1, n_observations))

        for index in range(0,n_observations):

            diffusive_flux_independent[:,index] = - mixture_mass_density[:,index] * np.dot(self.get_fickian_diffusion_coefficients_mass_mass[:,:,index], grad_species_mass_fractions[0:-1,index])

        diffusive_flux_last_species = - np.sum(diffusive_flux_independent, axis=0)

        diffusive_flux = np.vstack((diffusive_flux_independent, diffusive_flux_last_species[None,:]))

        return diffusive_flux

    # --------------------------------------------------------------------------

    def diffusive_flux_molar_volume(self, grad_species_molar_densities):
        """
        Computes the molar diffusive flux relative to a volume-averaged velocity, :math:`\\mathbf{J}_i^V`,
        from the generalized Fick's law:

        .. math::

            \\mathbf{J}_i^V = - \\mathbf{D}^V \\nabla \\mathbf{c}_i

        where :math:`\\mathbf{D}^V` are the molar Fickian diffusion coefficients expressed in a volume-averaged velocity reference frame.

        :param grad_species_molar_densities:
            vector ``numpy.ndarray`` specifying the species molar densities gradients, :math:`\\nabla \\mathbf{c}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.

        :return:
            - **diffusive_flux** - vector ``numpy.ndarray`` of molar diffusive fluxes relative to a volume-averaged velocity, :math:`\\mathbf{J}_i^V`, in :math:`[mole/(m^2s)]`. It has size ``(n_species,n_observations)``.
        """

        if self.get_fickian_diffusion_coefficients_molar_volume is None:
            raise ValueError("Parameter `fickian_diffusion_coefficients_molar_volume` has to be set in the `Diffusion` class object.")

        return diffusive_flux

    # --------------------------------------------------------------------------
