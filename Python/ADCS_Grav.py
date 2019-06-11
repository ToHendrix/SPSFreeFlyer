# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 16:39:44 2019

@author: J.J.M. Frijns
"""
# =============================================================================
# This program contains the calculation of the disturbance torque resulting 
# from the gravity field of the Earth
# =============================================================================

# =============================================================================
# Importing used packages, systems and files
# =============================================================================
from Database import magn_sph_data, E_to_P, mu
import numpy as np
from math import factorial as fact
from scipy.misc import derivative as der

# =============================================================================
# Defining arrays and variables
# =============================================================================


# =============================================================================
# Calculating gravity field potential
# =============================================================================

#n = 2
#m = 3
#delta = 3.1
#P_n = 1/((-2)**n*fact(n))*der(np.sin(x), x0=delta, dx=1e-6)
#P_nm = (1-np.(sin(delta))**2)**(m/2.)*der(np.sin(x), x0=delta, dx=1e-6)

M_ext_grav = 3