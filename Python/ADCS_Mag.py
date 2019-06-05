# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:32:59 2019

@author: J.J.M. Frijns
"""
# =============================================================================
# This program contains the calculation of the disturbance torque resulting 
# from the magnetic field of the Earth
# =============================================================================

# =============================================================================
# Importing used packages, systems and files
# =============================================================================
from Database import magn_data, coor_data
import numpy as np

# =============================================================================
# Defining arrays and variables
# =============================================================================
M_res_dip = np.array([[0.2],                    #Residual Dipole Moment [A m^2]
                    [0.2],
                    [0.2]])
M_ext_magn = np.zeros((43201,3))

# =============================================================================
# Initialising the arrays used
# =============================================================================

# Cropping the magnetic field array
magn_data = np.delete(magn_data, np.s_[0:6], axis=1)
magn_data = np.delete(magn_data, np.s_[-3:], axis=1)
for i in range(len(magn_data)):
    M_ext_magn[i] = np.cross(np.transpose(M_res_dip), np.array([magn_data[i]]))
