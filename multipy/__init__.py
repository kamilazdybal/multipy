"""multipy: Python library for multicomponent mass transfer"""

__author__ = "Kamila Zdybal"
__copyright__ = "Copyright (c) 2021, Kamila Zdybal"
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = ["Kamila Zdybal"]
__email__ = ["kamilazdybal@gmail.com"]
__status__ = "Production"

from .check import Check
from .composition import Composition
from .diffusion import Diffusion
from .flux import Flux
from .transform import Transform
from .velocity import Velocity
from .multicomponent_effects import MulticomponentEffects

from .plot import plot_species_fractions
from .plot import plot_1d_species_velocities
from .plot import plot_1d_averaged_velocities
from .plot import plot_1d_diffusive_flux
