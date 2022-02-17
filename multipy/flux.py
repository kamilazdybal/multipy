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
####    Class: Flux
####
################################################################################
################################################################################

class Flux:
    """
    Supports computing and storing fluxes. This class assumes that the species velocities, :math:`\\mathbf{u}_i`, are known.

    Diffusive fluxes:

    - mass diffusive flux relative to a mass-averaged velocity, :math:`\mathbf{j}_i`
    - mass diffusive flux relative to a molar-averaged velocity, :math:`\mathbf{j}_i^u`
    - molar diffusive flux relative to a mass-averaged velocity, :math:`\mathbf{J}_i^v`
    - molar diffusive flux relative to a molar-averaged velocity, :math:`\mathbf{J}_i`

    :param species_velocities:
        vector ``numpy.ndarray`` specifying the species velocities :math:`\mathbf{u}_i` in :math:`[m/s]`. It should be of size ``(n_species,n_observations)``.

    **Getters:**

    - **get_species_velocities**
    - **get_diffusive_molar_molar** (is set to ``None`` at class init)
    - **get_diffusive_molar_mass** (is set to ``None`` at class init)
    - **get_diffusive_mass_molar** (is set to ``None`` at class init)
    - **get_diffusive_mass_mass** (is set to ``None`` at class init)

    **Setters:**

    - **set_species_velocities**
    - **set_diffusive_molar_molar** (is set to ``None`` at class init)
    - **set_diffusive_molar_mass** (is set to ``None`` at class init)
    - **set_diffusive_mass_molar** (is set to ``None`` at class init)
    - **set_diffusive_mass_mass** (is set to ``None`` at class init)
    """

    # --------------------------------------------------------------------------

    def __init__(self, species_velocities):

        if not isinstance(species_velocities, np.ndarray):
            raise ValueError("Parameter `species_velocities` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_velocities)
        except:
            raise ValueError("Parameter `species_velocities` has to be a matrix.")

        if n_species < 2:
            raise ValueError("Parameter `species_velocities` has to have at least two species.")

        self.__species_velocities = species_velocities
        self.__velocity = multipy.Velocity(self.get_species_velocities)
        self.__diffusive_molar_molar = None
        self.__diffusive_molar_mass = None
        self.__diffusive_mass_molar = None
        self.__diffusive_mass_mass = None

    @property
    def get_species_velocities(self):
        return self.__species_velocities

    @property
    def get_diffusive_molar_molar(self):
        return self.__diffusive_molar_molar

    @property
    def get_diffusive_molar_mass(self):
        return self.__diffusive_molar_mass

    @property
    def get_diffusive_mass_molar(self):
        return self.__diffusive_mass_molar

    @property
    def get_diffusive_mass_mass(self):
        return self.__diffusive_mass_mass

    @get_species_velocities.setter
    def set_species_velocities(self, new_species_velocities):

        if new_species_velocities is not None:
            if not isinstance(new_species_velocities, np.ndarray):
                raise ValueError("Parameter `species_velocities` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_species_velocities)
            except:
                raise ValueError("Parameter `species_velocities` has to be a matrix.")

        self.__species_velocities = new_species_velocities

    @get_diffusive_molar_molar.setter
    def set_diffusive_molar_molar(self, new_diffusive_molar_molar):

        if new_diffusive_molar_molar is not None:
            if not isinstance(new_diffusive_molar_molar, np.ndarray):
                raise ValueError("Parameter `diffusive_molar_molar` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_diffusive_molar_molar)
            except:
                raise ValueError("Parameter `diffusive_molar_molar` has to be a matrix.")

        self.__diffusive_molar_molar = new_diffusive_molar_molar

    @get_diffusive_molar_mass.setter
    def set_diffusive_molar_mass(self, new_diffusive_molar_mass):

        if new_diffusive_molar_mass is not None:
            if not isinstance(new_diffusive_molar_mass, np.ndarray):
                raise ValueError("Parameter `diffusive_molar_mass` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_diffusive_molar_mass)
            except:
                raise ValueError("Parameter `diffusive_molar_mass` has to be a matrix.")

        self.__diffusive_molar_mass = new_diffusive_molar_mass

    @get_diffusive_mass_molar.setter
    def set_diffusive_mass_molar(self, new_diffusive_mass_molar):

        if new_diffusive_mass_molar is not None:
            if not isinstance(new_diffusive_mass_molar, np.ndarray):
                raise ValueError("Parameter `diffusive_mass_molar` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_diffusive_mass_molar)
            except:
                raise ValueError("Parameter `diffusive_mass_molar` has to be a matrix.")

        self.__diffusive_mass_molar = new_diffusive_mass_molar

    @get_diffusive_mass_mass.setter
    def set_diffusive_mass_mass(self, new_diffusive_mass_mass):

        if new_diffusive_mass_mass is not None:
            if not isinstance(new_diffusive_mass_mass, np.ndarray):
                raise ValueError("Parameter `diffusive_mass_mass` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_diffusive_mass_mass)
            except:
                raise ValueError("Parameter `diffusive_mass_mass` has to be a matrix.")

        self.__diffusive_mass_mass = new_diffusive_mass_mass

    # --------------------------------------------------------------------------

    def plot_diffusive_flux(self, species_names=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the computed diffusive fluxes.

        **Example:**

        .. image:: ../images/stefan-tube-diffusive-flux-molar-diff-molar-avg.svg
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

        if filename is not None:

            path = False

            if filename[0:2] == '..':
                __filename = filename[2::]
                path = True
            else:
                __filename = filename

            __base = __filename.split('.')[0]
            __extension = __filename.split('.')[1]

            if path:
                __filename = '..' + __base
            else:
                __filename = __base

        if self.get_diffusive_molar_molar is not None:
            if filename is not None:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_molar_molar, flux='molar', velocity='molar', species_names=species_names, colors=colors, figsize=figsize, filename=__filename + '-molar-diff-molar-avg.' + __extension)
            else:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_molar_molar, flux='molar', velocity='molar', species_names=species_names, colors=colors, figsize=figsize, filename=None)

        if self.get_diffusive_molar_mass is not None:
            if filename is not None:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_molar_mass, flux='molar', velocity='mass', species_names=species_names, colors=colors, figsize=figsize, filename=__filename + '-molar-diff-mass-avg.' + __extension)
            else:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_molar_mass, flux='molar', velocity='mass', species_names=species_names, colors=colors, figsize=figsize, filename=None)

        if self.get_diffusive_mass_molar is not None:
            if filename is not None:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_mass_molar, flux='mass', velocity='molar', species_names=species_names, colors=colors, figsize=figsize, filename=__filename + '-mass-diff-molar-avg.' + __extension)
            else:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_mass_molar, flux='mass', velocity='molar', species_names=species_names, colors=colors, figsize=figsize, filename=None)

        if self.get_diffusive_mass_mass is not None:
            if filename is not None:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_mass_mass, flux='mass', velocity='mass', species_names=species_names, colors=colors, figsize=figsize, filename=__filename + '-mass-diff-mass-avg.' + __extension)
            else:
                plt = multipy.plot.plot_1d_diffusive_flux(self.get_diffusive_mass_mass, flux='mass', velocity='mass', species_names=species_names, colors=colors, figsize=figsize, filename=None)

    # --------------------------------------------------------------------------

    def diffusive_molar_molar(self, species_mole_fractions, species_molar_densities):
        """
        Computes the molar diffusive flux relative to a molar-averaged velocity:

        .. math::

            \mathbf{J}_i = c_i \mathbf{u}_i + c_i \mathbf{u}

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying the species mole fractions, :math:`X_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param species_molar_densities:
            scalar ``numpy.ndarray`` specifying the molar densities of species, :math:`c_i`, in :math:`[mole/m^3]`. It should be of size ``(n_species,n_observations)``.

        :return:
            - **diffusive_flux** - vector ``numpy.ndarray`` of molar diffusive fluxes relative to a molar-averaged velocity :math:`\mathbf{J}_i` in :math:`[mole/(m^2s)]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_molar_densities, np.ndarray):
            raise ValueError("Parameter `species_molar_densities` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_molar_densities)
        except:
            raise ValueError("Parameter `species_molar_densities` has to be a matrix.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_molar_densities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_molar_densities` have different number of species `n_species`.")

        (n_species, n_observations) = np.shape(self.get_species_velocities)

        if n_observations != n_observations_1:
            raise ValueError("Parameters `species_mole_fractions`, `species_molar_densities` and `species_velocities` have different number of observations `n_observations`.")

        if n_species != n_species_1:
            raise ValueError("Parameters `species_mole_fractions`, `species_molar_densities` and `species_velocities` have different number of species `n_species`.")

        molar_averaged_velocity = self.__velocity.molar_averaged(species_mole_fractions)

        diffusive_flux = np.multiply(species_molar_densities, self.get_species_velocities) - np.multiply(species_molar_densities, molar_averaged_velocity)
        self.__diffusive_molar_molar = diffusive_flux

        return diffusive_flux

    # --------------------------------------------------------------------------

    def diffusive_molar_mass(self, species_mass_fractions, species_molar_densities):
        """
        Computes the molar diffusive flux relative to a mass-averaged velocity:

        .. math::

            \mathbf{J}_i^v = c_i \mathbf{u}_i + c_i \mathbf{v}

        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying the species mass fractions, :math:`Y_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param species_molar_densities:
            scalar ``numpy.ndarray`` specifying the species molar densities :math:`c_i` in :math:`[mole/m^3]`. It should be of size ``(n_species,n_observations)``.

        :return:
            - **diffusive_flux** - vector ``numpy.ndarray`` of molar diffusive fluxes relative to a mass-averaged velocity :math:`\mathbf{J}_i^v` in :math:`[mole/(m^2s)]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if not isinstance(species_molar_densities, np.ndarray):
            raise ValueError("Parameter `species_molar_densities` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_molar_densities)
        except:
            raise ValueError("Parameter `species_molar_densities` has to be a matrix.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_molar_densities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_molar_densities` have different number of species `n_species`.")

        (n_species, n_observations) = np.shape(self.get_species_velocities)

        if n_observations != n_observations_1:
            raise ValueError("Parameters `species_mass_fractions`, `species_molar_densities` and `species_velocities` have different number of observations `n_observations`.")

        if n_species != n_species_1:
            raise ValueError("Parameters `species_mass_fractions`, `species_molar_densities` and `species_velocities` have different number of species `n_species`.")

        mass_averaged_velocity = self.__velocity.mass_averaged(species_mass_fractions)

        diffusive_flux = np.multiply(species_molar_densities, self.get_species_velocities) - np.multiply(species_molar_densities, mass_averaged_velocity)
        self.__diffusive_molar_mass = diffusive_flux

        return diffusive_flux

    # --------------------------------------------------------------------------

    def diffusive_mass_molar(self, species_mole_fractions, species_mass_densities):
        """
        Computes the mass diffusive flux relative to a molar-averaged velocity:

        .. math::

            \mathbf{j}_i^u = \\rho_i \mathbf{u}_i + \\rho_i \mathbf{u}

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying the species mole fractions :math:`X_i` in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param species_mass_densities:
            scalar ``numpy.ndarray`` specifying the species mass densities :math:`\mathbf{\\rho}_i` in :math:`[kg/m^3]`. It should be of size ``(n_species,n_observations)``.

        :return:
            - **diffusive_flux** - vector ``numpy.ndarray`` of mass diffusive fluxes relative to a molar-averaged velocity :math:`\mathbf{j}_i^u` in :math:`[kg/(m^2s)]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(species_mass_densities, np.ndarray):
            raise ValueError("Parameter `species_mass_densities` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_mass_densities)
        except:
            raise ValueError("Parameter `species_mass_densities` has to be a matrix.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_densities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mole_fractions` and `species_mass_densities` have different number of species `n_species`.")

        (n_species, n_observations) = np.shape(self.get_species_velocities)

        if n_observations != n_observations_1:
            raise ValueError("Parameters `species_mole_fractions`, `species_mass_densities` and `species_velocities` have different number of observations `n_observations`.")

        if n_species != n_species_1:
            raise ValueError("Parameters `species_mole_fractions`, `species_mass_densities` and `species_velocities` have different number of species `n_species`.")

        molar_averaged_velocity = self.__velocity.molar_averaged(species_mole_fractions)

        diffusive_flux = np.multiply(species_mass_densities, self.get_species_velocities) - np.multiply(species_mass_densities, molar_averaged_velocity)
        self.__diffusive_mass_molar = diffusive_flux

        return diffusive_flux

    # --------------------------------------------------------------------------

    def diffusive_mass_mass(self, species_mass_fractions, species_mass_densities):
        """
        Computes the mass diffusive flux relative to a mass-averaged velocity:

        .. math::

            \mathbf{j}_i = \\rho_i \mathbf{u}_i + \\rho_i \mathbf{v}

        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying the species mass fractions :math:`Y_i` in :math:`[-]`. It should be of size ``(n_species, n_observations)``.
        :param species_mass_densities:
            scalar ``numpy.ndarray`` specifying the species mass densities :math:`\mathbf{\\rho}_i` in :math:`[kg/m^3]`. It should be of size ``(n_species, n_observations)``.

        :return:
            - **diffusive_flux** - vector ``numpy.ndarray`` of mass diffusive fluxes relative to a mass-averaged velocity :math:`\mathbf{j}_i` in :math:`[kg/(m^2s)]`. It has size ``(n_species, n_observations)``.
        """
        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if not isinstance(species_mass_densities, np.ndarray):
            raise ValueError("Parameter `species_mass_densities` has to be of type `numpy.ndarray`.")

        try:
            (n_species_2, n_observations_2) = np.shape(species_mass_densities)
        except:
            raise ValueError("Parameter `species_mass_densities` has to be a matrix.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_mass_densities` have different number of observations `n_observations`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_mass_fractions` and `species_mass_densities` have different number of species `n_species`.")

        (n_species, n_observations) = np.shape(self.get_species_velocities)

        if n_observations != n_observations_1:
            raise ValueError("Parameters `species_mass_fractions`, `species_mass_densities` and `species_velocities` have different number of observations `n_observations`.")

        if n_species != n_species_1:
            raise ValueError("Parameters `species_mass_fractions`, `species_mass_densities` and `species_velocities` have different number of species `n_species`.")

        mass_averaged_velocity = self.__velocity.mass_averaged(species_mass_fractions)

        diffusive_flux = np.multiply(species_mass_densities, self.get_species_velocities) - np.multiply(species_mass_densities, mass_averaged_velocity)
        self.__diffusive_mass_mass = diffusive_flux

        return diffusive_flux

    # --------------------------------------------------------------------------
