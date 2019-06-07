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
from Database import magn_sph_data, E_to_P
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
magn_sph_data = np.delete(magn_sph_data, np.s_[0:6], axis=1)
magn_sph_data = np.delete(magn_sph_data, np.s_[-3:], axis=1)
magn_sph_data[:,1] = -magn_sph_data[:,1]
magn_E_data = magn_sph_data

magn_P_data = E_to_P(magn_E_data)

for i in range(len(magn_P_data)):
    M_ext_magn[i] = np.cross(np.transpose(M_res_dip), [magn_P_data[i]])