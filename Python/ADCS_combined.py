# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:42:45 2019

@author: J.J.M. Frijns
"""

# =============================================================================
# This program contains the calculation of the combined disturbance torques
# resulting from the magnetic field of the Earth
# =============================================================================

# =============================================================================
# Importing used packages, systems and files
# =============================================================================
from ADCS_Grav import M_ext_grav
from ADCS_Mag import M_ext_magn
from ADCS_Aerodynamics import M_ext_aero
from Database import coor_elli_data
import numpy as np

# =============================================================================
# Defining arrays and variables
# =============================================================================
tau_build = np.zeros((43201,3))
tau_orb = np.zeros((431,4))
# =============================================================================
# Calculating the mission torque profile
# =============================================================================
M_ext = M_ext_grav + M_ext_magn + M_ext_aero
M_ext_x = np.delete(M_ext, np.s_[-2:], axis=1)
M_ext_y0 = np.delete(M_ext, np.s_[0], axis=1)
M_ext_y = np.delete(M_ext_y0, np.s_[-1], axis=1)
M_ext_z = np.delete(M_ext, np.s_[:2], axis=1)
M_ext_max = np.array([max(M_ext_x), max(M_ext_y), max(M_ext_z)])
M_ext_avg = np.array([np.sum(M_ext[:,0])/len(M_ext), \
                      np.sum(M_ext[:,1])/len(M_ext), \
                      np.sum(M_ext[:,2])/len(M_ext)])

ds = (coor_elli_data[1,0]-coor_elli_data[0,0])*86400.

for i in range(len(M_ext)):
    tau_build[i] =  M_ext[i]*ds
for i in range(int(len(tau_build)/100)-1):
    tau_orb[i,0], tau_orb[i,1], tau_orb[i,2], tau_orb[i,3] = \
    sum(tau_build[100*i:100*(i+1),0]) + sum(tau_build[100*i:100*(i+1),1]) + \
    sum(tau_build[100*i:100*(i+1),2]), sum(tau_build[100*i:100*(i+1),0]), \
    sum(tau_build[100*i:100*(i+1),1]), sum(tau_build[100*i:100*(i+1),2])       #momentum build-up per orbit (Total, x, y, z) [Nms]                          

tau_orb_max = np.array([np.amax(np.abs(tau_orb[:,1])), np.amax(np.abs(tau_orb[:,2])), \
                        np.amax(np.abs(tau_orb[:,3]))])
tau_max = np.amax(np.abs(tau_orb[:,0]))
