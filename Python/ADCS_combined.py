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
#from Database import
import numpy as np

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

#ds = 