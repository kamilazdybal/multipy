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
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------

def plot_species_fractions(species_fractions, fraction='mole', species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):

    __fraction = ['mole', 'mass', 'volume']

    if fraction not in __fraction:
        raise ValueError("Parameter `fraction` can only be 'molar', 'mass' or 'volume'.")

    try:
        (n_species, n_observations) = np.shape(species_fractions)
    except:
        raise ValueError("Parameter `species_fractions` has to be a matrix.")

    if species_names is not None:
        if not isinstance(species_names, list):
            raise ValueError("Parameter `species_names` has to be of type `list`.")

    if custom_coordinates is not None:
        if not isinstance(custom_coordinates, dict):
            raise ValueError("Parameter `custom_coordinates` has to be of type `dict`.")

    if colors is not None:
        if not isinstance(colors, list):
            raise ValueError("Parameter `colors` has to be of type `list`.")

    if not isinstance(figsize, tuple):
        raise ValueError("Parameter `figsize` has to be of type `tuple`.")

    if filename is not None:
        if not isinstance(filename, str):
            raise ValueError("Parameter `filename` has to be of type `str`.")

    if custom_coordinates is not None:
        x_label = next(iter(custom_coordinates))
        x = custom_coordinates[x_label]
    else:
        x_label = '$\eta$ [-]'
        x = np.linspace(0,1,n_observations)

    markevery=int(n_observations/10)

    fig = plt.figure(figsize=figsize)

    for i in range(0,n_species):
        if colors is not None:
            plt.plot(x, species_fractions[i,:], 'o-', color=colors[i], lw=1, markevery=markevery)
        else:
            plt.plot(x, species_fractions[i,:], 'o-', lw=1, markevery=markevery)

    plt.grid(alpha=0.3)
    plt.xlabel(x_label, fontsize=20)

    if fraction == 'mole':
        y_label = '$X_i$ [-]'
    elif fraction == 'mass':
        y_label = '$Y_i$ [-]'
    elif fraction == 'volume':
        y_label = '$V_i$ [-]'

    plt.ylabel(y_label, fontsize=20)
    plt.ylim([-0.05,1.05])
    plt.subplots_adjust(wspace=0, hspace=0.4)

    if species_names is not None:
        if n_species > 3:
            plt.legend(species_names, bbox_to_anchor=(1.05, 1))
        else:
            plt.legend(species_names)

    if filename is not None: plt.savefig(filename, dpi = 200, bbox_inches='tight')

    return plt

# ------------------------------------------------------------------------------

def plot_species_fractions_gradients(species_fractions_gradients, fraction='mole', species_names=None, custom_coordinates=None, colors=None, figsize=(10,5), filename=None):

    __fraction = ['mole', 'mass', 'volume']

    if fraction not in __fraction:
        raise ValueError("Parameter `fraction` can only be 'molar', 'mass' or 'volume'.")

    try:
        (n_species, n_observations) = np.shape(species_fractions_gradients)
    except:
        raise ValueError("Parameter `species_fractions_gradients` has to be a matrix.")

    if species_names is not None:
        if not isinstance(species_names, list):
            raise ValueError("Parameter `species_names` has to be of type `list`.")

    if custom_coordinates is not None:
        if not isinstance(custom_coordinates, dict):
            raise ValueError("Parameter `custom_coordinates` has to be of type `dict`.")

    if colors is not None:
        if not isinstance(colors, list):
            raise ValueError("Parameter `colors` has to be of type `list`.")

    if not isinstance(figsize, tuple):
        raise ValueError("Parameter `figsize` has to be of type `tuple`.")

    if filename is not None:
        if not isinstance(filename, str):
            raise ValueError("Parameter `filename` has to be of type `str`.")

    if custom_coordinates is not None:
        x_label = next(iter(custom_coordinates))
        x = custom_coordinates[x_label]
    else:
        x_label = '$\eta$ [-]'
        x = np.linspace(0,1,n_observations)

    markevery=int(n_observations/10)

    fig = plt.figure(figsize=figsize)

    for i in range(0,n_species):
        if colors is not None:
            plt.plot(x, species_fractions_gradients[i,:], 'o-', color=colors[i], lw=1, markevery=markevery)
        else:
            plt.plot(x, species_fractions_gradients[i,:], 'o-', lw=1, markevery=markevery)

    plt.grid(alpha=0.3)
    plt.xlabel('$\eta$ [-]', fontsize=20)

    if fraction == 'mole':
        y_label = r'$\nabla X_i$ [-]'
    elif fraction == 'mass':
        y_label = r'$\nabla Y_i$ [-]'
    elif fraction == 'volume':
        y_label = r'$\nabla V_i$ [-]'

    plt.ylabel(x_label, fontsize=20)
    # plt.ylim([-0.05,1.05])
    plt.subplots_adjust(wspace=0, hspace=0.4)

    if species_names is not None:
        if n_species > 3:
            plt.legend(species_names, bbox_to_anchor=(1.05, 1))
        else:
            plt.legend(species_names)

    if filename is not None: plt.savefig(filename, dpi = 200, bbox_inches='tight')

    return plt

# ------------------------------------------------------------------------------

def plot_1d_species_velocities(species_velocities, species_names=None, colors=None, figsize=(10,5), filename=None):

    try:
        (n_species, n_observations) = np.shape(species_velocities)
    except:
        raise ValueError("Parameter `species_velocities` has to be a matrix.")

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
            raise ValueError("Parameter `filename` has to be of type `str`.")

    x = np.linspace(0,1,n_observations)

    markevery=int(n_observations/10)

    fig = plt.figure(figsize=figsize)

    for i in range(0,n_species):
        if colors is not None:
            plt.plot(x, species_velocities[i,:], 'o-', color=colors[i], lw=1, markevery=markevery)
        else:
            plt.plot(x, species_velocities[i,:], 'o-', lw=1, markevery=markevery)

    plt.grid(alpha=0.3)
    plt.xlabel('$\eta$ [-]', fontsize=20)
    y_label = '$\mathbf{u}_i$ $[m/s]$'
    plt.ylabel(y_label, fontsize=20)
    plt.subplots_adjust(wspace=0, hspace=0.4)

    if species_names is not None:
        if n_species > 3:
            plt.legend(species_names, bbox_to_anchor=(1.05, 1))
        else:
            plt.legend(species_names)

    if filename is not None: plt.savefig(filename, dpi = 200, bbox_inches='tight')

    return plt

# ------------------------------------------------------------------------------

def plot_1d_averaged_velocities(velocity, colors=None, figsize=(10,5), filename=None):

    molar_averaged_velocity = velocity.get_molar_averaged
    mass_averaged_velocity = velocity.get_mass_averaged
    volume_averaged_velocity = velocity.get_volume_averaged

    if molar_averaged_velocity is not None: (_,n_observations) = np.shape(molar_averaged_velocity)
    if mass_averaged_velocity is not None: (_,n_observations) = np.shape(mass_averaged_velocity)
    if volume_averaged_velocity is not None: (_,n_observations) = np.shape(volume_averaged_velocity)

    try:
        x = np.linspace(0,1,n_observations)
    except:
        raise ValueError("No velocities to plot.")

    markevery=int(n_observations/10)

    fig = plt.figure(figsize=figsize)

    if colors is not None:
        if molar_averaged_velocity is not None: plt.plot(x, molar_averaged_velocity.ravel(), 'o-', color=colors[0], lw=1, markevery=markevery, label='Molar-averaged, $\mathbf{u}$')
        if mass_averaged_velocity is not None: plt.plot(x, mass_averaged_velocity.ravel(), 'o-', color=colors[1], lw=1, markevery=markevery, label='Mass-averaged, $\mathbf{v}$')
        if volume_averaged_velocity is not None: plt.plot(x, volume_averaged_velocity.ravel(), 'o-', color=colors[2], lw=1, markevery=markevery, label='Volume-averaged, $\mathbf{u}^V$')
    else:
        if molar_averaged_velocity is not None: plt.plot(x, molar_averaged_velocity.ravel(), 'o-', lw=1, markevery=markevery, label='Molar-averaged, $\mathbf{u}$')
        if mass_averaged_velocity is not None: plt.plot(x, mass_averaged_velocity.ravel(), 'o-', lw=1, markevery=markevery, label='Mass-averaged, $\mathbf{v}$')
        if volume_averaged_velocity is not None: plt.plot(x, volume_averaged_velocity.ravel(), 'o-', lw=1, markevery=markevery, label='Volume-averaged, $\mathbf{u}^V$')

    plt.grid(alpha=0.3)
    plt.xlabel('$\eta$ [-]', fontsize=20)
    plt.ylabel('Avg. velocity $[m/s]$', fontsize=20)
    plt.subplots_adjust(wspace=0, hspace=0.4)
    plt.legend()

    if filename is not None: plt.savefig(filename, dpi = 200, bbox_inches='tight')

    return plt

# ------------------------------------------------------------------------------

def plot_1d_diffusive_flux(diffusive_flux, flux='molar', velocity='molar', species_names=None, colors=None, figsize=(10,5), filename=None):

    __flux = ['molar', 'mass']
    __velocity = ['molar', 'mass']

    if not isinstance(diffusive_flux, np.ndarray):
        raise ValueError("Parameter `diffusive_flux` has to be of type `numpy.ndarray`.")

    try:
        (n_species, n_observations) = np.shape(diffusive_flux)
    except:
        raise ValueError("Parameter `diffusive_flux` has to be a matrix.")

    if not isinstance(flux, str):
        raise ValueError("Parameter `flux` has to be of type `str`.")

    if not isinstance(velocity, str):
        raise ValueError("Parameter `velocity` has to be of type `str`.")

    if flux not in __flux:
        raise ValueError("Parameter `flux` can only be 'molar' or 'mass'.")

    if velocity not in __velocity:
        raise ValueError("Parameter `velocity` can only be 'molar' or 'mass'.")

    if species_names is not None:
        if not isinstance(species_names, list):
            raise ValueError("Parameter `species_names` has to be of type `list`.")
        n_species_2 = len(species_names)
        if n_species != n_species_2:
            raise ValueError("Parameters `diffusive_flux` and `species_names` have different number of species `n_species`.")

    if not isinstance(figsize, tuple):
        raise ValueError("Parameter `figsize` has to be of type `tuple`.")

    if filename is not None:
        if not isinstance(filename, str):
            raise ValueError("Parameter `filename` has to be of type `str`.")

    x = np.linspace(0,1,n_observations)

    markevery=int(n_observations/10)

    fig = plt.figure(figsize=figsize)

    for i in range(0,n_species):
        if colors is not None:
            plt.plot(x, diffusive_flux[i,:], 'o-', color=colors[i], lw=1, markevery=markevery)
        else:
            plt.plot(x, diffusive_flux[i,:], 'o-', lw=1, markevery=markevery)

    plt.grid(alpha=0.3)
    plt.xlabel('$\eta$ [-]', fontsize=20)

    if flux == 'mass' and velocity == 'mass':
        y_label = '$\mathbf{j}_i$ $[kg/(m^2s)]$'
    elif flux == 'molar' and velocity == 'molar':
        y_label = '$\mathbf{J}_i$ $[mole/(m^2s)]$'
    elif flux == 'mass' and velocity == 'molar':
        y_label = '$\mathbf{j}_i^u$ $[kg/(m^2s)]$'
    elif flux == 'molar' and velocity == 'mass':
        y_label = '$\mathbf{J}_i^v$ $[mole/(m^2s)]$'

    plt.ylabel(y_label, fontsize=20)
    plt.subplots_adjust(wspace=0, hspace=0.4)

    if species_names is not None:
        if n_species > 3:
            plt.legend(species_names, bbox_to_anchor=(1.05, 1))
        else:
            plt.legend(species_names)

    if filename is not None: plt.savefig(filename, dpi = 200, bbox_inches='tight')

    return plt

# ------------------------------------------------------------------------------
