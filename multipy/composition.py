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
####    Class: Composition
####
################################################################################
################################################################################

class Composition:
    """
    Supports computing and storing composition-related quantities.

    **Getters:**

    - **get_species_mole_fractions** (is set to ``None`` at class init)
    - **get_species_mass_fractions** (is set to ``None`` at class init)
    - **get_species_volume_fractions** (is set to ``None`` at class init)
    - **get_grad_species_mole_fractions** (is set to ``None`` at class init)
    - **get_grad_species_mass_fractions** (is set to ``None`` at class init)
    - **get_grad_species_volume_fractions** (is set to ``None`` at class init)
    - **get_species_molar_densities** (is set to ``None`` at class init)
    - **get_species_mass_densities** (is set to ``None`` at class init)
    - **get_mixture_molar_density** (is set to ``None`` at class init)
    - **get_mixture_molar_volume** (is set to ``None`` at class init)
    - **get_mixture_molar_mass** (is set to ``None`` at class init)
    - **get_mixture_mass_density** (is set to ``None`` at class init)

    **Setters:**

    - **set_species_mole_fractions** setter for ``get_species_mole_fractions``
    - **set_species_mass_fractions** setter for ``get_species_mass_fractions``
    - **set_species_volume_fractions** setter for ``get_species_volume_fractions``
    - **set_grad_species_mole_fractions** setter for ``get_grad_species_mole_fractions``
    - **set_grad_species_mass_fractions** setter for ``get_grad_species_mass_fractions``
    - **set_grad_species_volume_fractions** setter for ``get_grad_species_volume_fractions``
    """

    # --------------------------------------------------------------------------

    def __init__(self):

        self.__species_mole_fractions = None
        self.__species_mass_fractions = None
        self.__species_volume_fractions = None

        self.__grad_species_mole_fractions = None
        self.__grad_species_mass_fractions = None
        self.__grad_species_volume_fractions = None

        self.__species_molar_densities = None
        self.__species_mass_densities = None

        self.__mixture_molar_density = None
        self.__mixture_molar_volume = None
        self.__mixture_mass_density = None
        self.__mixture_molar_mass = None

    @property
    def get_species_mole_fractions(self):
        return self.__species_mole_fractions

    @property
    def get_species_mass_fractions(self):
        return self.__species_mass_fractions

    @property
    def get_species_volume_fractions(self):
        return self.__species_volume_fractions

    @property
    def get_grad_species_mole_fractions(self):
        return self.__grad_species_mole_fractions

    @property
    def get_grad_species_mass_fractions(self):
        return self.__grad_species_mass_fractions

    @property
    def get_grad_species_volume_fractions(self):
        return self.__grad_species_volume_fractions

    @property
    def get_species_molar_densities(self):
        return self.__species_molar_densities

    @property
    def get_species_mass_densities(self):
        return self.__species_mass_densities





    @get_species_mole_fractions.setter
    def set_species_mole_fractions(self, new_species_mole_fractions):

        if new_species_mole_fractions is not None:
            if not isinstance(new_species_mole_fractions, np.ndarray):
                raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_species_mole_fractions)
            except:
                raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        self.__species_mole_fractions = new_species_mole_fractions

    @get_species_mass_fractions.setter
    def set_species_mass_fractions(self, new_species_mass_fractions):

        if new_species_mass_fractions is not None:
            if not isinstance(new_species_mass_fractions, np.ndarray):
                raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_species_mass_fractions)
            except:
                raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        self.__species_mass_fractions = new_species_mass_fractions

    @get_species_volume_fractions.setter
    def set_species_volume_fractions(self, new_species_volume_fractions):

        if new_species_volume_fractions is not None:
            if not isinstance(new_species_volume_fractions, np.ndarray):
                raise ValueError("Parameter `species_volume_fractions` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_species_volume_fractions)
            except:
                raise ValueError("Parameter `species_volume_fractions` has to be a matrix.")

        self.__species_volume_fractions = new_species_volume_fractions

    @get_grad_species_mole_fractions.setter
    def set_grad_species_mole_fractions(self, new_grad_species_mole_fractions):

        if new_grad_species_mole_fractions is not None:
            if not isinstance(new_grad_species_mole_fractions, np.ndarray):
                raise ValueError("Parameter `grad_species_mole_fractions` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_grad_species_mole_fractions)
            except:
                raise ValueError("Parameter `grad_species_mole_fractions` has to be a matrix.")

        self.__grad_species_mole_fractions = new_grad_species_mole_fractions

    @get_grad_species_mass_fractions.setter
    def set_grad_species_mass_fractions(self, new_grad_species_mass_fractions):

        if new_grad_species_mass_fractions is not None:
            if not isinstance(new_grad_species_mass_fractions, np.ndarray):
                raise ValueError("Parameter `grad_species_mass_fractions` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_grad_species_mass_fractions)
            except:
                raise ValueError("Parameter `grad_species_mass_fractions` has to be a matrix.")

        self.__grad_species_mass_fractions = new_grad_species_mass_fractions

    @get_grad_species_volume_fractions.setter
    def set_grad_species_volume_fractions(self, new_grad_species_volume_fractions):

        if new_grad_species_volume_fractions is not None:
            if not isinstance(new_grad_species_volume_fractions, np.ndarray):
                raise ValueError("Parameter `grad_species_volume_fractions` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_grad_species_volume_fractions)
            except:
                raise ValueError("Parameter `grad_species_volume_fractions` has to be a matrix.")

        self.__grad_species_volume_fractions = new_grad_species_volume_fractions

    @get_species_molar_densities.setter
    def set_species_molar_densities(self, new_species_molar_densities):

        if new_species_molar_densities is not None:
            if not isinstance(new_species_molar_densities, np.ndarray):
                raise ValueError("Parameter `species_molar_densities` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_species_molar_densities)
            except:
                raise ValueError("Parameter `species_molar_densities` has to be a matrix.")

        self.__species_molar_densities = new_species_molar_densities

    @get_species_mass_densities.setter
    def set_species_mass_densities(self, new_species_mass_densities):

        if new_species_mass_densities is not None:
            if not isinstance(new_species_mass_densities, np.ndarray):
                raise ValueError("Parameter `species_mass_densities` has to be of type `numpy.ndarray`.")

            try:
                (n_species, n_observations) = np.shape(new_species_mass_densities)
            except:
                raise ValueError("Parameter `species_mass_densities` has to be a matrix.")

        self.__species_mass_densities = new_species_mass_densities

    # --------------------------------------------------------------------------

    def plot_species_mole_fractions(self, species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the computed species mole fractions, :math:`\\mathbf{X}_i`.

        **Example:**

        .. image:: ../images/species-mole-fractions.svg
          :width: 400

        :param species_names: (optional)
            ``list`` of ``str`` specifying the species names.
        :param custom_coordinates: (optional)
            ``dict`` specifying the custom coordinates on the :math:`x` axis. It should be of the format ``{ label : array }``
            where ``label`` is a string that will be plotted as an :math:`x` axis label and ``array`` is a vector
            defining the custom grid. If not specified, a generic uniform grid between 0 and 1 will be used.
        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each species. Example: ``colors=['#C7254E', '#BBBBBB', '#008CBA']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.

        :return:
            - **plot_handle** - ``matplotlib.pyplot`` plot handle.
        """

        if self.get_species_mole_fractions is not None:
            plot_handle = multipy.plot.plot_species_fractions(self.get_species_mole_fractions, fraction='mole', species_names=species_names, custom_coordinates=custom_coordinates, colors=colors, figsize=figsize, filename=filename)

            return plot_handle

    # --------------------------------------------------------------------------

    def plot_species_mass_fractions(self, species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the computed species mass fractions, :math:`\\mathbf{Y}_i`.

        **Example:**

        .. image:: ../images/species-mass-fractions.svg
          :width: 400

        :param species_names: (optional)
            ``list`` of ``str`` specifying the species names.
        :param custom_coordinates: (optional)
            ``dict`` specifying the custom coordinates on the :math:`x` axis. It should be of the format ``{ label : array }``
            where ``label`` is a string that will be plotted as an :math:`x` axis label and ``array`` is a vector
            defining the custom grid. If not specified, a generic uniform grid between 0 and 1 will be used.
        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each species. Example: ``colors=['#C7254E', '#BBBBBB', '#008CBA']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.

        :return:
            - **plot_handle** - ``matplotlib.pyplot`` plot handle.
        """

        if self.get_species_mass_fractions is not None:
            plot_handle = multipy.plot.plot_species_fractions(self.get_species_mass_fractions, fraction='mass', species_names=species_names, custom_coordinates=custom_coordinates, colors=colors, figsize=figsize, filename=filename)

            return plot_handle

    # --------------------------------------------------------------------------

    def plot_species_volume_fractions(self, species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the computed species volume fractions, :math:`\\mathbf{V}_i`.

        **Example:**

        .. image:: ../images/species-volume-fractions.svg
          :width: 400

        :param species_names: (optional)
            ``list`` of ``str`` specifying the species names.
        :param custom_coordinates: (optional)
            ``dict`` specifying the custom coordinates on the :math:`x` axis. It should be of the format ``{ label : array }``
            where ``label`` is a string that will be plotted as an :math:`x` axis label and ``array`` is a vector
            defining the custom grid. If not specified, a generic uniform grid between 0 and 1 will be used.
        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each species. Example: ``colors=['#C7254E', '#BBBBBB', '#008CBA']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.

        :return:
            - **plot_handle** - ``matplotlib.pyplot`` plot handle.
        """

        if self.get_species_volume_fractions is not None:
            plot_handle = multipy.plot.plot_species_fractions(self.get_species_volume_fractions, fraction='volume', species_names=species_names, custom_coordinates=custom_coordinates, colors=colors, figsize=figsize, filename=filename)

            return plot_handle

    # --------------------------------------------------------------------------

    def plot_grad_species_mole_fractions(self, species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the computed species mole fraction gradients, :math:`\\nabla \\mathbf{X}_i`.

        **Example:**

        .. image:: ../images/species-mole-fraction-gradients.svg
          :width: 400

        :param species_names: (optional)
            ``list`` of ``str`` specifying the species names.
        :param custom_coordinates: (optional)
            ``dict`` specifying the custom coordinates on the :math:`x` axis. It should be of the format ``{ label : array }``
            where ``label`` is a string that will be plotted as an :math:`x` axis label and ``array`` is a vector
            defining the custom grid. If not specified, a generic uniform grid between 0 and 1 will be used.
        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each species. Example: ``colors=['#C7254E', '#BBBBBB', '#008CBA']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.

        :return:
            - **plot_handle** - ``matplotlib.pyplot`` plot handle.
        """

        if self.get_grad_species_mole_fractions is not None:
            plot_handle = multipy.plot.plot_species_fractions_gradients(self.get_grad_species_mole_fractions, fraction='mole', species_names=species_names, custom_coordinates=custom_coordinates, colors=colors, figsize=figsize, filename=filename)

            return plot_handle

    # --------------------------------------------------------------------------

    def plot_grad_species_mass_fractions(self, species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the computed species mass fraction gradients, :math:`\\nabla \\mathbf{Y}_i`.

        **Example:**

        .. image:: ../images/species-mass-fraction-gradients.svg
          :width: 400

        :param species_names: (optional)
            ``list`` of ``str`` specifying the species names.
        :param custom_coordinates: (optional)
            ``dict`` specifying the custom coordinates on the :math:`x` axis. It should be of the format ``{ label : array }``
            where ``label`` is a string that will be plotted as an :math:`x` axis label and ``array`` is a vector
            defining the custom grid. If not specified, a generic uniform grid between 0 and 1 will be used.
        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each species. Example: ``colors=['#C7254E', '#BBBBBB', '#008CBA']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.

        :return:
            - **plot_handle** - ``matplotlib.pyplot`` plot handle.
        """

        if self.get_grad_species_mass_fractions is not None:
            plot_handle = multipy.plot.plot_species_fractions_gradients(self.get_grad_species_mass_fractions, fraction='mass', species_names=species_names, custom_coordinates=custom_coordinates, colors=colors, figsize=figsize, filename=filename)

            return plot_handle

    # --------------------------------------------------------------------------

    def plot_grad_species_volume_fractions(self, species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):
        """
        Plots the computed species volume fraction gradients, :math:`\\nabla \\mathbf{V}_i`.

        **Example:**

        .. image:: ../images/species-volume-fraction-gradients.svg
          :width: 400

        :param species_names: (optional)
            ``list`` of ``str`` specifying the species names.
        :param custom_coordinates: (optional)
            ``dict`` specifying the custom coordinates on the :math:`x` axis. It should be of the format ``{ label : array }``
            where ``label`` is a string that will be plotted as an :math:`x` axis label and ``array`` is a vector
            defining the custom grid. If not specified, a generic uniform grid between 0 and 1 will be used.
        :param colors: (optional)
            ``list`` of ``str`` specifying the plotting colors for each species. Example: ``colors=['#C7254E', '#BBBBBB', '#008CBA']``.
        :param figsize: (optional)
            ``tuple`` specifying the figure size.
        :param filename: (optional)
            ``str`` specifying the filename. If set to ``None``, plot will not be saved to a file.

        :return:
            - **plot_handle** - ``matplotlib.pyplot`` plot handle.
        """

        if self.get_grad_species_volume_fractions is not None:
            plot_handle = multipy.plot.plot_species_fractions_gradients(self.get_grad_species_volume_fractions, fraction='volume', species_names=species_names, custom_coordinates=custom_coordinates, colors=colors, figsize=figsize, filename=filename)

            return plot_handle

    # --------------------------------------------------------------------------

    def mixture_molar_density(self, T, p):
        """
        Computes the mixture molar density:

        .. math::

            c = p/RT

        under the assumption of the ideal gas law, where the mixture molar density,
        :math:`c`, is constant at constant temperature, :math:`T`, and pressure, :math:`p`.

        :param T:
            ``int`` or ``float`` specifying the temperature, :math:`T`, in :math:`[K]`.
        :param p:
            ``int`` or ``float`` specifying the pressure, :math:`p`, in :math:`[Pa]`.

        :return:
            - **mixture_molar_density** - mixture molar density, :math:`c`, in :math:`[mole/m^3]`.
        """

        if not isinstance(T, int) and not isinstance(T, float):
                raise ValueError("Parameter `T` has to be of type `int` or `float`.")

        if not isinstance(p, int) and not isinstance(p, float):
            raise ValueError("Parameter `p` has to be of type `int` or `float`.")

        if T < 0:
            raise ValueError("Parameter `T` cannot be negative.")

        if p < 0:
            raise ValueError("Parameter `p` cannot be negative.")
            
        mixture_molar_density = p/(gas_constant*T)

        return mixture_molar_density

    # --------------------------------------------------------------------------

    def mixture_mass_density(self, species_mass_densities):
        """
        Computes the mixture mass density:

        .. math::

            \\rho = \\sum_{i=1}^{n} \\rho_i

        :param species_mass_densities:
            scalar ``numpy.ndarray`` specifying **all** species mass densities, :math:`\\pmb{\\rho}_i`, in :math:`[kg/m^3]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.

        :return:
            - **mixture_mass_density** - scalar ``numpy.ndarray`` of mixture mass density, :math:`\\pmb{\\rho}`, in :math:`[kg/m^3]`. It has size ``(1,n_observations)``.
        """

        if not isinstance(species_mass_densities, np.ndarray):
            raise ValueError("Parameter `species_mass_densities` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_mass_densities)
        except:
            raise ValueError("Parameter `species_mass_densities` has to be a matrix.")

        if n_species < 2:
            raise ValueError("Parameter `species_mass_densities` should contain all species. Only one species found.")

        mixture_mass_density = np.sum(species_mass_densities, axis=0)[np.newaxis,:]

        return mixture_mass_density

    # --------------------------------------------------------------------------

    def mixture_molar_mass(self, species_fractions, basis, species_molar_masses):
        """
        Computes the mixture molar mass at each observation, either from the species mole fractions
        (when ``basis='molar'``):

        .. math::

            \\mathbf{M} = \\sum_{i=1}^{n} X_i M_i

        or from species mass fractions (when ``basis='mass'``):

        .. math::

            \\mathbf{M} = \\Big( \\sum_{i=1}^{n} Y_i M_i \\Big)^{-1}

        :param species_fractions:
            scalar ``numpy.ndarray`` specifying **all** species mole (or mass) fractions, :math:`X_i` (or :math:`Y_i`), in :math:`[-]`.
            It should be of size ``(n_species,n_observations)`` where ``n_species`` is at least 2.
        :param basis:
            ``str`` specifying whether the molar or mass equation should be used. Can be ``'molar'`` or ``'mass'``.
        :param species_molar_masses:
            scalar ``numpy.ndarray`` specifying the species molar masses, :math:`\\mathbf{M}_i`, in :math:`[kg/mole]`.
            It should be of size ``(n_species,1)``.

        :return:
            - **mixture_molar_mass** - scalar ``numpy.ndarray`` of mixture molar masses, :math:`\\pmb{M}`, in :math:`[kg/mole]`. It has size ``(1,n_observations)``.
        """

        __basis = ['molar', 'mass']

        if not isinstance(species_fractions, np.ndarray):
            raise ValueError("Parameter `species_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations_1) = np.shape(species_fractions)
        except:
            raise ValueError("Parameter `species_fractions` has to be a matrix.")

        if n_species_1 < 2:
            raise ValueError("Parameter `species_fractions` has to have at least two species.")

        if not isinstance(species_molar_masses, np.ndarray):
            raise ValueError("Parameter `species_molar_masses` has to be of type `numpy.ndarray`.")

        if basis not in __basis:
            raise ValueError("Parameter `basis` has to be 'molar' or 'mass'.")

        try:
            (n_species_2, n_dim) = np.shape(species_molar_masses)
        except:
            raise ValueError("Parameter `species_molar_masses` has to be a matrix.")

        if n_dim != 1:
            raise ValueError("Parameter `species_molar_masses` has to be of size ``(n_species,1)``.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_fractions` and `species_molar_masses` have different number of species, ``n_species``.")

        if np.any(species_molar_masses==0):
            raise ValueError("Parameter `species_molar_masses` has entries equal to zero.")

        mixture_molar_mass = np.zeros((1,n_observations_1))

        if basis == 'molar':
            for i in range(0, n_observations_1):
                mixture_molar_mass[0,i] = np.sum(np.multiply(species_fractions[:,i][:,None], species_molar_masses))
        elif basis == 'mass':
            for i in range(0, n_observations_1):
                mixture_molar_mass [0,i] = 1.0 / np.sum(np.divide(species_fractions[:,i][:,None], species_molar_masses))

        return mixture_molar_mass

    # --------------------------------------------------------------------------

    def mixture_molar_volume(self, T, p):
        """
        Computes the mixture molar volume, :math:`\\bar{V}`, from:

        .. math::

            \\bar{V} = \\frac{1}{c}

        under the assumption of the ideal gas law, where the mixture molar density,
        :math:`c`, is constant at constant temperature, :math:`T`, and pressure, :math:`p`:

        .. math::

            c = p/RT

        :param T: (optional)
            ``int`` or ``float`` specifying the temperature, :math:`T`, in :math:`[K]`.
        :param p: (optional)
            ``int`` or ``float`` specifying the pressure, :math:`p`, in :math:`[Pa]`.

        :return:
            - **mixture_molar_volume** - mixture molar volume, :math:`\\bar{V}`, in :math:`[m^3/mole]`.
        """

        if not isinstance(T, int) and not isinstance(T, float):
                raise ValueError("Parameter `T` has to be of type `int` or `float`.")

        if not isinstance(p, int) and not isinstance(p, float):
            raise ValueError("Parameter `p` has to be of type `int` or `float`.")

        mixture_molar_density = self.mixture_molar_density(T, p)

        mixture_molar_volume = 1 / mixture_molar_density
        self.__mixture_molar_volume = mixture_molar_volume

        return mixture_molar_volume

    # --------------------------------------------------------------------------

    def species_mole_fractions(self, species_molar_densities, T, p):
        """
        Computes the species mole fractions:

        .. math::

            X_i = \\frac{c_i}{c}

        under the assumption of the ideal gas law, where the mixture molar density,
        :math:`c`, is constant at constant temperature, :math:`T`, and pressure, :math:`p`:

        .. math::

            c = p/RT

        :param species_molar_densities:
            scalar ``numpy.ndarray`` specifying the species molar densities, :math:`\\mathbf{c}_i`, in :math:`[mole/m^3]`.
            It should be of size ``(n_species,n_observations)``. Note that ``n_species`` can be equal to 1,
            if you are computing the mole fraction for only one of the species in the mixture.
        :param T:
            ``int`` or ``float`` specifying the temperature, :math:`T`, in :math:`[K]`.
        :param p:
            ``int`` or ``float`` specifying the pressure, :math:`p`, in :math:`[Pa]`.

        :return:
            - **species_mole_fractions** - scalar ``numpy.ndarray`` of species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_molar_densities, np.ndarray):
            raise ValueError("Parameter `species_molar_densities` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_molar_densities)
        except:
            raise ValueError("Parameter `species_molar_densities` has to be a matrix.")

        if not isinstance(T, int) and not isinstance(T, float):
                raise ValueError("Parameter `T` has to be of type `int` or `float`.")

        if not isinstance(p, int) and not isinstance(p, float):
            raise ValueError("Parameter `p` has to be of type `int` or `float`.")

        mixture_molar_density = self.mixture_molar_density(T, p)

        species_mole_fractions = species_molar_densities / mixture_molar_density
        self.__species_mole_fractions = species_mole_fractions

        return species_mole_fractions

    # --------------------------------------------------------------------------

    def species_mass_fractions(self, species_mass_densities, mixture_mass_density):
        """
        Computes the species mass fractions:

        .. math::

            Y_i = \\frac{\\rho_i}{\\rho}

        :param species_mass_densities:
            scalar ``numpy.ndarray`` specifying the species mass densities, :math:`\\pmb{\\rho}_i`, in :math:`[kg/m^3]`.
            It should be of size ``(n_species,n_observations)``. Note that ``n_species`` can be equal to 1,
            if you are computing the mass fraction for only one of the species in the mixture.
        :param mixture_mass_density:
            scalar ``numpy.ndarray`` specifying the mixture mass density, :math:`\\pmb{\\rho}`, in :math:`[kg/m^3]`.
            It should be of size ``(1,n_observations)``.

        :return:
            - **species_mass_fractions** - scalar ``numpy.ndarray`` of species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mass_densities, np.ndarray):
            raise ValueError("Parameter `species_mass_densities` has to be of type `numpy.ndarray`.")

        if not isinstance(mixture_mass_density, np.ndarray):
            raise ValueError("Parameter `mixture_mass_density` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations_1) = np.shape(species_mass_densities)
        except:
            raise ValueError("Parameter `species_mass_densities` has to be a matrix.")

        (n_dim, n_observations_2) = np.shape(mixture_mass_density)

        if n_dim != 1:
            raise ValueError("Parameter `mixture_mass_density` has to be of size `(1, n_observations)`.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_mass_densities` and `mixture_mass_density` have different number of observations, `n_observations`.")

        species_mass_fractions = np.divide(species_mass_densities, mixture_mass_density)
        self.__species_mass_fractions = species_mass_fractions

        return species_mass_fractions

    # --------------------------------------------------------------------------

    def species_volume_fractions(self, species_volume, mixture_volume):
        """
        Computes the species mass fractions:

        .. math::

            V_i = \\frac{v_i}{V}

        :param species_volume:
            scalar ``numpy.ndarray`` specifying the species volumes, :math:`v_i`, in :math:`[m^3]`.
            It should be of size ``(n_species,n_observations)``. Note that ``n_species`` can be equal to 1,
            if you are computing the volume fraction for only one of the species in the mixture.
        :param mixture_volume:
            scalar ``numpy.ndarray`` specifying the mixture volume, :math:`V`, in :math:`[m^3]`.
            It should be of size ``(1,n_observations)``.

        :return:
            - **species_volume_fractions** - scalar ``numpy.ndarray`` of species volume fractions, :math:`\\mathbf{V}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_volume, np.ndarray):
            raise ValueError("Parameter `species_volume` has to be of type `numpy.ndarray`.")

        if not isinstance(mixture_volume, np.ndarray):
            raise ValueError("Parameter `mixture_volume` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations_1) = np.shape(species_volume)
        except:
            raise ValueError("Parameter `species_volume` has to be a matrix.")

        (n_dim, n_observations_2) = np.shape(mixture_volume)

        if n_dim != 1:
            raise ValueError("Parameter `mixture_volume` has to be of size `(1, n_observations)`.")

        if n_observations_1 != n_observations_2:
            raise ValueError("Parameters `species_volume` and `mixture_volume` have different number of observations, `n_observations`.")

        species_volume_fractions = np.divide(species_volume, mixture_volume)
        self.__species_volume_fractions = species_volume_fractions

        return species_volume_fractions

    # --------------------------------------------------------------------------

    def species_molar_densities(self, species_mole_fractions, T, p):
        """
        Computes the species molar densities:

        .. math::

            c_i = c X_i

        under the assumption of the ideal gas law, where the mixture molar density,
        :math:`c`, is constant at constant temperature, :math:`T`, and pressure, :math:`p`:

        .. math::

            c = p/RT

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying the species mole fractions :math:`X_i` in :math:`[-]`.
            It should be of size ``(n_observations,n_species)``. Note that ``n_species`` can be equal to 1,
            if you are computing the molar density for only one of the species in the mixture.
        :param T:
            ``int`` or ``float`` specifying the temperature :math:`T` in :math:`[K]`.
        :param p:
            ``int`` or ``float`` specifying the pressure :math:`p` in :math:`[Pa]`.

        :return:
            - **species_molar_densities** - scalar ``numpy.ndarray`` of species molar densities :math:`c_i` in :math:`[mole/m^3]`. It has size ``(n_observations,n_species)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_observations, n_species) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(T, int) and not isinstance(T, float):
                raise ValueError("Parameter `T` has to be of type `int` or `float`.")

        if not isinstance(p, int) and not isinstance(p, float):
            raise ValueError("Parameter `p` has to be of type `int` or `float`.")

        mixture_molar_density = self.mixture_molar_density(T, p)

        species_molar_densities = species_mole_fractions * mixture_molar_density
        self.__species_molar_densities = species_molar_densities

        return species_molar_densities

    # --------------------------------------------------------------------------

    def species_mass_densities(self, species_molar_densities, species_molar_masses):
        """
        Computes the species mass densities:

        .. math::

            \mathbf{\\rho}_i = c_i M_i

        :param species_molar_densities:
            scalar ``numpy.ndarray`` specifying the species molar densities, :math:`\\mathbf{c}_i`, in :math:`[mole/m^3]`.
            It should be of size ``(n_species,n_observations)``. Note that ``n_species`` can be equal to 1,
            if you are computing the mass density for only one of the species in the mixture.
        :param species_molar_masses:
            scalar ``numpy.ndarray`` specifying the species molar masses, :math:`\\mathbf{M}_i`, in :math:`[kg/mole]`.
            It should be of size ``(n_species,1)``.

        :return:
            - **species_mass_densities** - scalar ``numpy.ndarray`` of species mass densities, :math:`\\pmb{\\mathbf{\\rho}}_i`, in :math:`[kg/m^3]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_molar_densities, np.ndarray):
            raise ValueError("Parameter `species_molar_densities` has to be of type `numpy.ndarray`.")

        if not isinstance(species_molar_masses, np.ndarray):
            raise ValueError("Parameter `species_molar_masses` has to be of type `numpy.ndarray`.")

        try:
            (n_species_1, n_observations) = np.shape(species_molar_densities)
        except:
            raise ValueError("Parameter `species_molar_densities` has to be a matrix.")

        try:
            (n_species_2, n_dim) = np.shape(species_molar_masses)
        except:
            raise ValueError("Parameter `species_molar_masses` has to be a matrix.")

        if n_dim != 1:
            raise ValueError("Parameter `species_molar_masses` has to be of size `(n_species,1)`.")

        if n_species_1 != n_species_2:
            raise ValueError("Parameters `species_molar_densities` and `species_molar_masses` have different number of species `n_species`.")

        if np.any(species_molar_masses==0):
            raise ValueError("Parameter `species_molar_masses` has entries equal to zero.")

        species_mass_densities = np.multiply(species_molar_densities, species_molar_masses)

        self.__species_mass_densities = species_mass_densities

        return species_mass_densities

    # --------------------------------------------------------------------------

    def grad_species_mole_fractions(self, species_mole_fractions, delta, edge_order=1):
        """
        Computes species mole fraction gradients, :math:`\\nabla \\mathbf{X}_i`,
        numerically from the species mole fractions, :math:`X_i`:

        .. math::

            \\nabla X_i = \\partial_d X_i

        assuming a uniform grid spacing in a spatial dimension :math:`d`.

        .. note::

            ``numpy.gradient`` is used to compute gradients. Gradients are
            computed using second order accurate central differences in the interior
            points and either first or second order accurate
            (forward or backward) differences at the boundaries.

        :param species_mole_fractions:
            scalar ``numpy.ndarray`` specifying the species mole fractions, :math:`\\mathbf{X}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param delta:
            ``int`` or ``float`` specifying the grid spacing in the dimension :math:`d`.
        :param edge_order:
            ``int`` specifying whether first or second order accurate differences are computed at the boundaries. It should be ``1`` or ``2``.

        :return:
            - **species_mole_fractions_gradients** - vector ``numpy.ndarray`` of species mole fractions gradients, :math:`\\nabla \\mathbf{X}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mole_fractions, np.ndarray):
            raise ValueError("Parameter `species_mole_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_mole_fractions)
        except:
            raise ValueError("Parameter `species_mole_fractions` has to be a matrix.")

        if not isinstance(delta, int) and not isinstance(delta, float):
            raise ValueError("Parameter `delta` has to be of type `int` or `float`.")

        if edge_order != 1 and edge_order != 2:
            raise ValueError("Parameter `edge_order` can only be `1` or `2`.")

        species_mole_fractions_gradients = np.zeros_like(species_mole_fractions)

        for i in range(0,n_species):

            d_species_mole_fractions = np.gradient(species_mole_fractions[i,:], edge_order=edge_order)
            species_mole_fractions_gradients[i,:] = d_species_mole_fractions / delta

        self.__grad_species_mole_fractions = species_mole_fractions_gradients

        return species_mole_fractions_gradients

    # --------------------------------------------------------------------------

    def grad_species_mass_fractions(self, species_mass_fractions, delta, edge_order=1):
        """
        Computes species mass fraction gradients, :math:`\\nabla \\mathbf{Y}_i`,
        numerically from the species mass fractions, :math:`Y_i`:

        .. math::

            \\nabla Y_i = \\partial_d Y_i

        assuming a uniform grid spacing in a spatial dimension :math:`d`.

        .. note::

            ``numpy.gradient`` is used to compute gradients. Gradients are
            computed using second order accurate central differences in the interior
            points and either first or second order accurate
            (forward or backward) differences at the boundaries.

        :param species_mass_fractions:
            scalar ``numpy.ndarray`` specifying the species mass fractions, :math:`\\mathbf{Y}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param delta:
            ``int`` or ``float`` specifying the grid spacing in the dimension :math:`d`.
        :param edge_order:
            ``int`` specifying whether first or second order accurate differences are computed at the boundaries. It should be ``1`` or ``2``.

        :return:
            - **species_mass_fractions_gradients** - vector ``numpy.ndarray`` of species mass fractions gradients, :math:`\\nabla \\mathbf{Y}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_mass_fractions, np.ndarray):
            raise ValueError("Parameter `species_mass_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_mass_fractions)
        except:
            raise ValueError("Parameter `species_mass_fractions` has to be a matrix.")

        if not isinstance(delta, int) and not isinstance(delta, float):
            raise ValueError("Parameter `delta` has to be of type `int` or `float`.")

        if edge_order != 1 and edge_order != 2:
            raise ValueError("Parameter `edge_order` can only be `1` or `2`.")

        species_mass_fractions_gradients = np.zeros_like(species_mass_fractions)

        for i in range(0,n_species):

            d_species_mass_fractions = np.gradient(species_mass_fractions[i,:], edge_order=edge_order)
            species_mass_fractions_gradients[i,:] = d_species_mass_fractions / delta

        self.__grad_species_mass_fractions = species_mass_fractions_gradients

        return species_mass_fractions_gradients

    # --------------------------------------------------------------------------

    def grad_species_volume_fractions(self, species_volume_fractions, delta, edge_order=1):
        """
        Computes species volume fraction gradients, :math:`\\nabla \\mathbf{V}_i`,
        numerically from the species volume fractions, :math:`V_i`:

        .. math::

            \\nabla V_i = \\partial_d V_i

        assuming a uniform grid spacing in a spatial dimension :math:`d`.

        .. note::

            ``numpy.gradient`` is used to compute gradients. Gradients are
            computed using second order accurate central differences in the interior
            points and either first or second order accurate
            (forward or backward) differences at the boundaries.

        :param species_volume_fractions:
            scalar ``numpy.ndarray`` specifying the species volume fractions, :math:`\\mathbf{V}_i`, in :math:`[-]`. It should be of size ``(n_species,n_observations)``.
        :param delta:
            ``int`` or ``float`` specifying the grid spacing in the dimension :math:`d`.
        :param edge_order:
            ``int`` specifying whether first or second order accurate differences are computed at the boundaries. It should be ``1`` or ``2``.

        :return:
            - **species_volume_fractions_gradients** - vector ``numpy.ndarray`` of species volume fractions gradients, :math:`\\nabla \\mathbf{V}_i`, in :math:`[-]`. It has size ``(n_species,n_observations)``.
        """

        if not isinstance(species_volume_fractions, np.ndarray):
            raise ValueError("Parameter `species_volume_fractions` has to be of type `numpy.ndarray`.")

        try:
            (n_species, n_observations) = np.shape(species_volume_fractions)
        except:
            raise ValueError("Parameter `species_volume_fractions` has to be a matrix.")

        if not isinstance(delta, int) and not isinstance(delta, float):
            raise ValueError("Parameter `delta` has to be of type `int` or `float`.")

        if edge_order != 1 and edge_order != 2:
            raise ValueError("Parameter `edge_order` can only be `1` or `2`.")

        species_volume_fractions_gradients = np.zeros_like(species_volume_fractions)

        for i in range(0,n_species):

            d_species_volume_fractions = np.gradient(species_volume_fractions[i,:], edge_order=edge_order)
            species_volume_fractions_gradients[i,:] = d_species_volume_fractions / delta

        self.__grad_species_volume_fractions = species_volume_fractions_gradients

        return species_volume_fractions_gradients
