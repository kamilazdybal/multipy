"""multipy: Python library for multicomponent mass transfer"""

__author__ = "James C. Sutherland, Kamila Zdybal"
__copyright__ = "Copyright (c) 2021, James C. Sutherland, Kamila Zdybal"
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
####    Class: Velocity
####
################################################################################
################################################################################

class Velocity:
    """
    Supports computing and storing velocities.

    Species convective velocities:

    - species velocities, :math:`\\mathbf{u}_i`

    Mixture-averaged convective velocities:

    - mass-averaged velocity, :math:`\\mathbf{v}`
    - molar-averaged velocity, :math:`\\mathbf{u}`
    - volume-averaged velocity, :math:`\\mathbf{u}^V`
    - arbitrarily-averaged velocity, :math:`\\mathbf{u}^a`

    :param species_velocities: (optional)
        vector ``numpy.ndarray`` specifying the species velocities :math:`\mathbf{u}_i`
        in :math:`[m/s]`. It should be of size ``(n_observations,n_species)``
        where ``n_species`` is at least 2.

    **Getters:**

    - **get_species_velocities**
    - **get_molar_averaged** (is set to ``None`` at class init)
    - **get_mass_averaged** (is set to ``None`` at class init)
    - **get_volume_averaged** (is set to ``None`` at class init)
    - **get_arbitrarily_averaged** (is set to ``None`` at class init)

    **Setters:**

    - **set_species_velocities** setter for ``get_species_velocities``
    - **set_molar_averaged** setter for ``get_molar_averaged``
    - **set_mass_averaged** setter for ``get_mass_averaged``
    - **set_volume_averaged** setter for ``get_volume_averaged``
    - **set_arbitrarily_averaged** setter for ``get_arbitrarily_averaged``
    """

    # --------------------------------------------------------------------------

    def __init__(self, species_velocities=None):

        if species_velocities is not None:
            if not isinstance(species_velocities, np.ndarray):
                raise ValueError("Parameter `species_velocities` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(species_velocities)
            except:
                raise ValueError("Parameter `species_velocities` has to be a matrix.")

            if n_species < 2:
                raise ValueError("Parameter `species_velocities` has to have at least two species.")

        self.__species_velocities = species_velocities
        self.__molar_averaged = None
        self.__mass_averaged = None
        self.__volume_averaged = None
        self.__arbitrarily_averaged = None

    @property
    def get_species_velocities(self):
        return self.__species_velocities

    @property
    def get_molar_averaged(self):
        return self.__molar_averaged

    @property
    def get_mass_averaged(self):
        return self.__mass_averaged

    @property
    def get_volume_averaged(self):
        return self.__volume_averaged

    @property
    def get_arbitrarily_averaged(self):
        return self.__arbitrarily_averaged

    @get_species_velocities.setter
    def set_species_velocities(self, new_species_velocities):

        if new_species_velocities is not None:
            if not isinstance(new_species_velocities, np.ndarray):
                raise ValueError("Parameter `species_velocities` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_species_velocities)
            except:
                raise ValueError("Parameter `species_velocities` has to be a matrix.")

            if n_species < 2:
                raise ValueError("Parameter `species_velocities` has to have at least two species.")

        self.__species_velocities = new_species_velocities

    @get_molar_averaged.setter
    def set_molar_averaged(self, new_molar_averaged):

        if new_molar_averaged is not None:
            if not isinstance(new_molar_averaged, np.ndarray):
                raise ValueError("Parameter `molar_averaged` has to be of type `numpy.ndarray`.")

            try:
                (n_dim, n_observations) = np.shape(new_molar_averaged)
            except:
                raise ValueError("Parameter `molar_averaged` has to have size `(1,n_observations)`.")

            if n_dim != 1:
                raise ValueError("Parameter `molar_averaged` has to have size `(1,n_observations)`.")

        self.__molar_averaged = new_molar_averaged

    @get_mass_averaged.setter
    def set_mass_averaged(self, new_mass_averaged):

        if new_mass_averaged is not None:
            if not isinstance(new_mass_averaged, np.ndarray):
                raise ValueError("Parameter `mass_averaged` has to be of type `numpy.ndarray`.")

            try:
                (n_dim, n_observations) = np.shape(new_mass_averaged)
            except:
                raise ValueError("Parameter `mass_averaged` has to have size `(1,n_observations)`.")

            if n_dim != 1:
                raise ValueError("Parameter `mass_averaged` has to have size `(1,n_observations)`.")

        self.__mass_averaged = new_mass_averaged

    @get_volume_averaged.setter
    def set_volume_averaged(self, new_volume_averaged):

        if new_volume_averaged is not None:
            if not isinstance(new_volume_averaged, np.ndarray):
                raise ValueError("Parameter `volume_averaged` has to be of type `numpy.ndarray`.")

            try:
                (n_dim, n_observations) = np.shape(new_volume_averaged)
            except:
                raise ValueError("Parameter `volume_averaged` has to have size `(1,n_observations)`.")

            if n_dim != 1:
                raise ValueError("Parameter `volume_averaged` has to have size `(1,n_observations)`.")

        self.__volume_averaged = new_volume_averaged

    @get_arbitrarily_averaged.setter
    def set_arbitrarily_averaged(self, new_arbitrarily_averaged):

        if new_arbitrarily_averaged is not None:
            if not isinstance(new_arbitrarily_averaged, np.ndarray):
                raise ValueError("Parameter `new_arbitrarily_averaged` has to be of type `numpy.ndarray`.")

            try:
                (n_dim, n_observations) = np.shape(new_arbitrarily_averaged)
            except:
                raise ValueError("Parameter `arbitrarily_averaged` has to have size `(1,n_observations)`.")

            if n_dim != 1:
                raise ValueError("Parameter `arbitrarily_averaged` has to have size `(1,n_observations)`.")

        self.__arbitrarily_averaged = new_arbitrarily_averaged

    # --------------------------------------------------------------------------

    def plot_species_velocities(self, species_names=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the species velocities.

        **Example:**

        .. image:: ../images/stefan-tube-species-velocities.svg
          :width: 400

        :param species_names: (optional)
            ``list`` of ``str`` specifying the species names.
        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each species. Example: ``colors=['#C7254E', '#BBBBBB', '#008CBA']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.
        """

        if species_names is not None:
            if not isinstance(species_names, list):
                raise ValueError("Parameter `species_names` has to be of type `list`.")

        if colors is not None:
            if not isinstance(colors, list):
                raise ValueError("Parameter `colors` has to be of type `list`.")

        if not isinstance(figsize, tuple):
            raise ValueError("Parameter `figsize` has to be of type `tuple`.")

        if filename is not None:
            if not isinstance(filename, str):
                raise ValueError("Parameter `filename` has to be of type `list`.")

        (n_species, n_observations) = np.shape(self.get_species_velocities)

        if n_species != len(species_names):
            raise ValueError("Parameter `species_names` has different number of species than `Velocity.get_species_velocities`.")

        if species_names is not None and colors is not None:
            if len(species_names) > len(colors):
                raise ValueError("Not enough colors specified for all the species.")

        multipy.plot.plot_1d_species_velocities(self.get_species_velocities, species_names=species_names, colors=colors, figsize=figsize, filename=filename)

    # --------------------------------------------------------------------------

    def plot_averaged_velocities(self, colors=None, figsize=(10,5), filename=None):
        """
        Plots the averaged velocities.

        **Example:**

        .. image:: ../images/stefan-tube-averaged-velocities.svg
          :width: 400

        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each averaged velocity. Example: ``colors=['#222222', '#858585']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.
        """

        if colors is not None:
            if not isinstance(colors, list):
                raise ValueError("Parameter `colors` has to be of type `list`.")

        if not isinstance(figsize, tuple):
            raise ValueError("Parameter `figsize` has to be of type `tuple`.")

        if filename is not None:
            if not isinstance(filename, str):
                raise ValueError("Parameter `filename` has to be of type `list`.")

        n_velocities = 0

        if self.get_molar_averaged is not None:
            n_velocities += 1
        if self.get_mass_averaged is not None:
            n_velocities += 1
        if self.get_volume_averaged is not None:
            n_velocities += 1

        if n_velocities > len(colors):
            raise ValueError("Not enough colors specified for all the averaged velocities.")

        multipy.plot.plot_1d_averaged_velocities(self, colors=colors, figsize=figsize, filename=filename)

    # --------------------------------------------------------------------------

    def species_velocities(self, total_flux, species_fractions, basis='molar', mixture_molar_density=None, mixture_mass_density=None):
        """
        Computes the species velocities, :math:`\\mathbf{u}_i`.

        If ``basis`` is set to ``'molar'``, species velocities are computed
        using the total molar fluxes, :math:`\\mathbf{N}_i`, species mole fractions, :math:`X_i`
        and the mixture molar density, :math:`c`:

        .. math::

            \\mathbf{u}_i = \\frac{\\mathbf{N}_i}{c X_i}

        If ``basis`` is set to ``'mass'``, species velocities are computed
        using the total mass fluxes, :math:`\\mathbf{n}_i`, species mass fractions, :math:`Y_i`
        and the mixture mass density, :math:`\\rho`:

        .. math::

            \\mathbf{u}_i = \\frac{\\mathbf{n}_i}{\\rho Y_i}

        :param total_flux:
            vector ``numpy.ndarray`` of total molar fluxes :math:`\\mathbf{N}_i`
            in :math:`[mole/(m^2s)]` or total mass fluxes :math:`\\mathbf{n}_i` in :math:`[kg/(m^2s)]`.
            It should be of size ``(n_species, n_observations)``.
        :param species_fractions:
            scalar ``numpy.ndarray`` specifying the species mole fractions :math:`X_i` in :math:`[-]` or species mass fractions :math:`Y_i` in :math:`[-]`. It should be of size ``(n_species, n_observations)``.
        :param basis: (optional)
            ``str`` specifying whether the molar or mass total flux equation should be used. Can be ``'molar'`` or ``'mass'``.
        :param mixture_molar_density: (optional)
            mixture molar density :math:`c` in :math:`[mole/m^3]`. Has to be specified if ``basis`` is set to ``molar``.
        :param mixture_mass_density: (optional)
            mixture mass density :math:`\\rho` in :math:`[kg/m^3]`. Has to be specified if ``basis`` is set to ``mass``.

        :return:
            - **species_velocities** - vector ``numpy.ndarray`` of species velocities :math:`\mathbf{u}_i` in :math:`[m/s]`. It has size ``(n_species, n_observations)``.
        """

        __basis = ['molar', 'mass']

        if not isinstance(total_flux, np.ndarray):
            raise ValueError("Parameter `total_flux` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(total_flux)
        except:
            raise ValueError("Parameter `total_flux` has to be a matrix.")

        if not isinstance(species_fractions, np.ndarray):
            raise ValueError("Parameter `species_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_fractions)
        except:
            raise ValueError("Parameter `species_fractions` has to be a matrix.")

        if basis not in __basis:
            raise ValueError("Parameter `basis` has to be 'molar' or 'mass'.")

        if basis == 'molar':
            if mixture_molar_density is None:
                raise ValueError("Parameter `mixture_molar_density` has to specified when basis is 'molar'.")

        if basis == 'mass':
            if mixture_mass_density is None:
                raise ValueError("Parameter `mixture_mass_density` has to specified when basis is 'mass'.")

        if basis == 'molar':
            species_velocities = np.divide(total_flux, np.multiply(mixture_molar_density, species_fractions))
        elif basis == 'mass':
            species_velocities = np.divide(total_flux, np.multiply(mixture_mass_density, species_fractions))

        self.__species_velocities = species_velocities

        return species_velocities

    # --------------------------------------------------------------------------

    def molar_averaged(self, species_mole_fractions):
        """
        Computes the molar-averaged velocity:

        .. math::

            \mathbf{u} = \\sum_{i=1}^{n} X_i \mathbf{u}_i

        where :math:`n` is the number of species.

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole fractions :math:`X_i` in :math:`[-]`. It should be of size ``(n_species,n_observations)`` where
            ``n_species`` is at least 2.

        :return:
            - **molar_averaged_velocity** - vector ``numpy.ndarray`` of molar-averaged velocity :math:`\mathbf{u}` in :math:`[m/s]`. It has size ``(1,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if n_species_1 < 2:
            raise ValueError("Parameter `species_mole_fractions` has to have at least two species.")

        (n_species_2, n_observations_2) = np.shape(self.get_species_velocities)

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_velocities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_velocities` have different number of species `n_species`.")

        idx = multipy.Check().sum_of_species_fractions(species_mole_fractions, tolerance=0.001)

        if len(idx) != 0:
            print('Warning: species mole fractions do not sum up to 1.0 within a 0.001 tolerance.')

        molar_averaged_velocity = np.sum(np.multiply(species_mole_fractions, self.get_species_velocities), axis=0)[None,:]
        self.__molar_averaged = molar_averaged_velocity

        return molar_averaged_velocity

    # --------------------------------------------------------------------------

    def mass_averaged(self, species_mass_fractions):
        """
        Computes the mass-averaged velocity:

        .. math::

            \mathbf{v} = \\sum_{i=1}^{n} Y_i \mathbf{u}_i

        where :math:`n` is the number of species.

        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mass fractions :math:`Y_i` in :math:`[-]`. It should be of size ``(n_species,n_observations)`` where
            ``n_species`` is at least 2.

        :return:
            - **mass_averaged_velocity** - vector ``numpy.ndarray`` of mass-averaged velocity in :math:`[m/s]`. It has size ``(1,n_observations)``.
        """

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_observations_1, n_species_1) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if n_species_1 < 2:
            raise ValueError("Parameter `species_mass_fractions` has to have at least two species.")

        (n_observations_2, n_species_2) = np.shape(self.get_species_velocities)

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_velocities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_velocities` have different number of species `n_species`.")

        idx = multipy.Check().sum_of_species_fractions(species_mass_fractions, tolerance=0.001)

        if len(idx) != 0:
            print('Warning: species mass fractions do not sum up to 1.0 within a 0.001 tolerance.')

        mass_averaged_velocity = np.sum(np.multiply(species_mass_fractions, self.get_species_velocities), axis=0)[None,:]
        self.__mass_averaged = mass_averaged_velocity

        return mass_averaged_velocity

    # --------------------------------------------------------------------------

    def volume_averaged(self, species_volume_fractions):
        """
        Computes the volume-averaged velocity:

        .. math::

            \mathbf{u}^V = \\sum_{i=1}^{n} V_i \mathbf{u}_i

        where :math:`n` is the number of species.

        :param species_volume_fractions:
            scalar ``numpy.ndarray`` specifying **all** species volume fractions :math:`V_i` in :math:`[-]`. It should be of size ``(n_species,n_observations)`` where
            ``n_species`` is at least 2.

        :return:
            - **volume_averaged_velocity** - vector ``numpy.ndarray`` of volume-averaged velocity :math:`\\mathbf{u}^V` in :math:`[m/s]`. It has size ``(1,n_observations)``.
        """

        if not isinstance(species_volume_fractions, np.ndarray):
            raise ValueError("Parameter `species_volume_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_volume_fractions)
        except:
            raise ValueError("Parameter `species_volume_fractions` has to be a matrix.")

        if n_species_1 < 2:
            raise ValueError("Parameter `species_volume_fractions` has to have at least two species.")

        (n_species_2, n_observations_2) = np.shape(self.get_species_velocities)

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_volume_fractions` and `species_velocities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_volume_fractions` and `species_velocities` have different number of species `n_species`.")

        idx = multipy.Check().sum_of_species_fractions(species_volume_fractions, tolerance=0.001)

        if len(idx) != 0:
            print('Warning: species volume fractions do not sum up to 1.0 within a 0.001 tolerance.')

        volume_averaged_velocity = np.sum(np.multiply(species_volume_fractions, self.get_species_velocities), axis=0)[None,:]
        self.__volume_averaged = volume_averaged_velocity

        return volume_averaged_velocity

    # --------------------------------------------------------------------------

    def arbitrarily_averaged(self, arbitrary_weighting_factors):
        """
        Computes the arbitrarily-averaged velocity:

        .. math::

            \mathbf{u}^a = \\sum_{i=1}^{n} a_i \mathbf{u}_i

        where :math:`n` is the number of species and :math:`a_i` are the arbitrary weighting factors, such that :math:`\\sum_{i=1}^{n} a_i = 1`.

        :param arbitrary_weighting_factors:
            scalar ``numpy.ndarray`` specifying arbitrary weighting factors, :math:`a_i` in :math:`[-]`, for **all** species. It should be of size ``(n_species,n_observations)`` where
            ``n_species`` is at least 2.

        :return:
            - **arbitrarily_averaged_velocity** - vector ``numpy.ndarray`` of arbitrarily-averaged velocity :math:`\mathbf{u}^a` in :math:`[m/s]`. It has size ``(1,n_observations)``.
        """

        if not isinstance(arbitrary_weighting_factors, np.ndarray):
            raise ValueError("Parameter `arbitrary_weighting_factors` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(arbitrary_weighting_factors)
        except:
            raise ValueError("Parameter `arbitrary_weighting_factors` has to be a matrix.")

        if n_species_1 < 2:
            raise ValueError("Parameter `arbitrary_weighting_factors` has to have at least two species.")

        (n_species_2, n_observations_2) = np.shape(self.get_species_velocities)

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `arbitrary_weighting_factors` and `species_velocities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `arbitrary_weighting_factors` and `species_velocities` have different number of species `n_species`.")

        idx = multipy.Check().sum_of_species_fractions(arbitrary_weighting_factors, tolerance=0.001)

        if len(idx) != 0:
            print('Warning: species volume fractions do not sum up to 1.0 within a 0.001 tolerance.')

        arbitrarily_averaged_velocity = np.sum(np.multiply(arbitrary_weighting_factors, self.get_species_velocities), axis=0)[None,:]
        self.__arbitrarily_averaged = arbitrarily_averaged_velocity

        return arbitrarily_averaged_velocity
