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
####    Class: MulticomponentEffects
####
################################################################################
################################################################################

class MulticomponentEffects:
    """
    Supports performing checks of the multicomponent effects.
    """

    # --------------------------------------------------------------------------

    def __init__(self):

        pass

    # --------------------------------------------------------------------------
