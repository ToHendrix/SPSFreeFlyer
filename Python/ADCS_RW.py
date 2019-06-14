# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:42:45 2019

@author: J.J.M. Frijns
"""

# =============================================================================
# This program contains the analysis of the different possible momentum wheels 
# used in the design of the SPS-2 FF
# =============================================================================

# =============================================================================
# Importing used packages, systems and files
# =============================================================================
from ADCS_combined import tau_orb
import numpy as np
# https://www.researchgate.net/profile/Abolfazl_Shirazi2/publication/276216468_Pyramidal_reaction_wheel_arrangement_optimization_of_satellite_attitude_control_subsystem_for_minimizing_power_consumption/links/566e4c6d08ae1a797e4061cc.pdf?origin=publication_list

# =============================================================================
# 
# =============================================================================


tau = 1.
kappa = np.deg2rad(32)
tau_y = 0.25*np.sin(kappa)