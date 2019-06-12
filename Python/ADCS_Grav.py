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
from Database import mu_earth, coor_P_data, m, r_out, h
import numpy as np

# =============================================================================
# Defining arrays and variables
# =============================================================================
M_ext_grav = np.zeros((43201,3))

# =============================================================================
# Calculating MMOIs
# =============================================================================
I_xx = 1./12.*m*(3*r_out**2 + h**2)
I_yy = I_xx
I_zz = 0.5*m*r_out**2

# =============================================================================
# Calculating gravity field potential
# =============================================================================
#Source: https://link.springer.com/content/pdf/bbm%3A978-3-642-25749-0%2F1.pdf
for i in range(len(coor_P_data)):
    r_p = coor_P_data[i,:3]/(np.sqrt((coor_P_data[i,0])**2+(coor_P_data[i,1])**2\
        +(coor_P_data[i,2])**2))
    J = np.array([[I_xx, 0, 0],
                  [0, I_yy, 0],
                  [0, 0, I_zz]])
    M_ext_grav[i] = list(3*mu_earth/np.sqrt((coor_P_data[i,0])**2+(coor_P_data[i,1])**2\
            +(coor_P_data[i,2])**2)**3 * np.cross(r_p, np.dot(J, r_p)))  