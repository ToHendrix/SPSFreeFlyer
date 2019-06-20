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
from Database import magn_sph_data, E_to_P1, P1_to_P, z_P, Sun_E_data, \
coor_P1_data , coor_E_data
import numpy as np

# =============================================================================
# Defining arrays and variables
# =============================================================================
M_res_dip = np.array([[0.2],                                                   #Residual Dipole Moment [A m^2]
                    [0.2],
                    [0.2]])
M_ext_magn = np.zeros((43201,3))                                               #Empty array for the magnetic field moment

# =============================================================================
# Initialising the arrays used
# =============================================================================
magn_sph_data = np.delete(magn_sph_data, np.s_[0:6], axis=1)                   #Cropping the magnetic field array to get spherical coordinates
magn_sph_data = np.delete(magn_sph_data, np.s_[-3:], axis=1)                   #Cropping the magnetic field array to get spherical coordinates 
magn_sph_data[:,1] = -magn_sph_data[:,1]                                       #Changing the direction of the soutward magnetic field component to be pointing north 
magn_sph_data[:,0] = -magn_sph_data[:,0] 
magn_E_data = magn_sph_data*10**-4                                             #Writing the magnetic field array to the correct named variable and unit [T]

magn_P1_data = E_to_P1(magn_E_data, Sun_E_data, coor_E_data)                   #Tranforming the magnetic field from the E- to P-frame         
magn_P_data = P1_to_P(magn_P1_data, z_P, coor_P1_data)

for i in range(len(magn_P_data)):
    M_ext_magn[i] = np.cross(magn_P_data[i], np.transpose(M_res_dip))          #Calculating magnetic field disturbance torque at every instance 