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
from ADCS_combined import tau_orb, tau_orb_max
import numpy as np
# https://www.researchgate.net/profile/Abolfazl_Shirazi2/publication/276216468_Pyramidal_reaction_wheel_arrangement_optimization_of_satellite_attitude_control_subsystem_for_minimizing_power_consumption/links/566e4c6d08ae1a797e4061cc.pdf?origin=publication_list

# =============================================================================
# Defining arrays and variables
# =============================================================================
tau_RW = np.zeros((431,3))
RW_req_tau = np.zeros((1,3))
RW_cap_mom = np.array([[1.1],
                   [0.03],
                   [0.011076],
                   [0.2],
                   [0.5],
                   [1],
                   [1],
                   [2.00E-02],
                   [1]])
RW_orb = np.zeros((9,2))

# =============================================================================
# Calculating the amount of orbits possible per wheel type
# =============================================================================
kappa = np.deg2rad(32)

RW_req_tau = [0.5*np.cos(kappa)*tau_orb_max[0], \
          0.5*np.cos(kappa)*tau_orb_max[1], 0.25*np.sin(kappa)*tau_orb_max[2]]
for i in range(len(tau_orb)):
    tau_RW[i] = [0.5*np.cos(kappa)*tau_orb[i,1], \
          0.5*np.cos(kappa)*tau_orb[i,2], 0.25*np.sin(kappa)*tau_orb[i,3]]
    
for i in range(len(RW_cap_mom)):
    RW_orb[i,0] = np.max(np.abs([np.sum(tau_RW[:,0])/RW_cap_mom[i], np.sum(tau_RW[:,1])/RW_cap_mom[i], \
          np.sum(tau_RW[:,2])/RW_cap_mom[i]]))
    RW_orb[i,1] = 431./RW_orb[i,0] 