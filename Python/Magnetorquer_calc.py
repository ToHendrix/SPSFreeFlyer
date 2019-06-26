# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:18:09 2019

@author: tom-h_000
"""

from ADCS_Mag import magn_P_data
from ADCS_combined import tau_orb
import numpy as np

magn_P_crop = np.transpose(np.array([0.5*(magn_P_data[:,0] + magn_P_data[:,1]), magn_P_data[:,2]]))

tau_orb_crop = np.transpose(np.array([tau_orb[:,1],  tau_orb[:,2], tau_orb[:,3]]))

# =============================================================================
# Function to calculate the angle of the magnetic field vector wrt the 
# P-frame. The data is given in x,y,z components in the magn_data list
# =============================================================================
def angle_magn(magnetic_data):
    #angle wrt the x axis
    theta_B = np.arctan(magnetic_data[1]/magnetic_data[0])
    #angle wrt the y-axis, this should be pi/2 rad together with the theta
    omega_B = np.arctan(magnetic_data[0]/magnetic_data[1])
    #angle wrt the z-axis
    phi_B = np.arctan(np.sqrt(magnetic_data[0]**2 + magnetic_data[1]**2)/ magnetic_data[2])
    


    return theta_B, omega_B, phi_B

# =============================================================================
# Calculate the torque generated from the magnetorquers, it is scaled by the 
# dipole vector which depends on the angle wrt the magnetic field
# =============================================================================
def get_Js(dipole, count):
<<<<<<< HEAD
    orbit_torque = dipole * abs(magn_P_data[count])
=======
    orbit_torque = dipole * abs(magn_P_data[count]) *600
#    print("orbit torque \n", orbit_torque)
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
    return orbit_torque

    
