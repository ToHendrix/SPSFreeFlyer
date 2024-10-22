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
from ADCS_combined import tau_orb, tau_orb_max, M_ext_avg
import numpy as np
import matplotlib.pyplot as plt
# https://www.researchgate.net/profile/Abolfazl_Shirazi2/publication/276216468_Pyramidal_reaction_wheel_arrangement_optimization_of_satellite_attitude_control_subsystem_for_minimizing_power_consumption/links/566e4c6d08ae1a797e4061cc.pdf?origin=publication_list

# =============================================================================
# Defining arrays and variables
# =============================================================================
tau_bld_RW = np.zeros((431,3))

RW_cap_mom = np.array([[1.1],
                   [0.03],
                   [0.011076],
                   [0.2],
                   [0.5],
                   [1],
                   [1],
                   [2.00E-02],
                   [1]])
RW_tau_orb = np.zeros((431,1))
RW_num_orb = np.zeros((9,2))

# =============================================================================
# Calculating the amount of orbits possible per wheel type
# =============================================================================
kappa = np.deg2rad(32)                                                         #Angle of reaction wheels w.r.t. (x,y) plane in the P-frame (rad)     

RW_req_tau = [0.5*tau_orb_max[0]*np.cos(kappa), \
          0.5*tau_orb_max[1]*np.cos(kappa), 0.25*tau_orb_max[2]*np.sin(kappa)] #Required torque per momentum wheel in each axis (x, y, z) [Nm]
RW_avg_tau = [0.5*M_ext_avg[0]*np.cos(kappa), \
          0.5*M_ext_avg[1]*np.cos(kappa), 0.25*M_ext_avg[2]*np.sin(kappa)]

for i in range(len(tau_orb)):
    tau_bld_RW[i] = [0.5*np.cos(kappa)*tau_orb[i,1], \
          0.5*np.cos(kappa)*tau_orb[i,2], 0.25*np.sin(kappa)*tau_orb[i,3]]     #Momentum build-up of a reaction wheel in the pyramid system per orbit (x, y, z) [Nm]
    RW_tau_orb[i] = sum(tau_bld_RW[i])                                         #Total momentum build-up of a reaction wheel in the pyramid system per orbit [Nm]

for i in range(len(RW_num_orb)):
    RW_num_orb[i] = [sum(np.abs(RW_tau_orb))/RW_cap_mom[i], \
              431/(sum(np.abs(RW_tau_orb))/RW_cap_mom[i])]                     #Amount of momentum dumps in 30 days, Amount of orbits before a momentum dump 
    

def calculate_system_charasteristics(tau_max):
    tau_max_wheel = tau_max                                                         #Max torque delivered by one momentum wheel. [Nm]
    tau_z = 4 * tau_max_wheel * np.sin(kappa) 
    tau_x = 2 * tau_max_wheel * np.cos(kappa)
    tau_y = 2 * tau_max_wheel * np.cos(kappa)
    
    return tau_x, tau_y, tau_z

def leave_eclipse(tau_x, tau_y, tau_z):
    initial_angle = [45 * (np.pi/180), 0 * (np.pi/180), 0* (np.pi/180)]
    
    angle = initial_angle
    omega = [0, 0, 0]
    count = 0
    
    MMOI = [45.868, 45.868, 83.298]
    
    omega_dot = [tau_x /  MMOI[0], tau_y / MMOI[1], tau_z / MMOI[2]]
    
    while (angle[0]  > initial_angle[0]/2) or (angle[1]  > initial_angle[1]/2) or (angle[2]  > initial_angle[2]/2):
        for i in range(3):
            omega[i] += omega_dot[i]
            
            angle[i] -= omega[i]
            
        count += 1
    
    return count *2

#turningtime = []
#x_values = []
#for i in range(1, 11):
#    j = i/100
#    print(j)
#    torques = calculate_system_charasteristics(j)
#    turningtime.append(leave_eclipse(torques[0], torques[1], torques[2]))
#    x_values.append(j)
#    plt.plot(x_values, turningtime)
#    plt.title("Time to regain sun pointing attitude after eclipse for different reaction wheels")
#    plt.ylabel("Time to change angle(min)")
#    plt.xlabel("Torque Cap of the wheel [Nm]")
#    plt.show()
    
torques = calculate_system_charasteristics(0.1)
print(leave_eclipse(torques[0], torques[1], torques[2]))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    