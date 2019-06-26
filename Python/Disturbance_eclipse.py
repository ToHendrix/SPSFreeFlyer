# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:53:23 2019

@author: tom-h_000
"""

import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
from ADCS_combined import tau_build
=======
from ADCS_combined import tau_orb
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
from Precession_influence import calculate_detumble_time

# =============================================================================
# Constants
# =============================================================================

t_ecl = 35                                                                     #eclipse time [min]

MMOI = [45.868,
       45.868,
       83.298]                                                                  #kg * m^2

mass = 250                                                                      #kg
# =============================================================================
# From matlab the values for b for the orbit will have to be imported to get
# the most accurate results possible
# =============================================================================

file = open('bees.txt', "r")
lines = file.readlines()                                                        #Read lines 
file.close()
b = []
for s in lines:                                 
    b.append(s.replace('i', '').replace('0-', '').replace('0+', '').replace(',', ' ').split())   

bees = np.asarray(b).astype(np.float)

def calculate_angular_velocity():
    #create empty matrix in which the omegas will be added, as it is assumed 
    #the sc is stable before eclipse the first entry will be null
    omega_mat = []
    omega_mat.append([0, 0, 0])
    attitude = ([0,0,0])
<<<<<<< HEAD
    omega_dot = [0,0,0]                                                         #initial acceleration, assumed to be 0
=======
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
    
    for minute in range(0, t_ecl):    
        # =============================================================================
        # Calculation of J, the inertia dyadic
        # =============================================================================
        
        #calculating the entries of the J matrix, as defined in the book\cite
        rho_1 = np.sqrt(MMOI[1]/mass + MMOI[0]/mass + MMOI[2]/mass)
        rho_3 = np.sqrt(abs(MMOI[1]/mass - rho_1**2))
        rho_2 = np.sqrt(abs(MMOI[0]/mass - rho_3**2))
        
        J_mat = [[MMOI[0],-rho_1*rho_2*mass, -rho_1*rho_3*mass],
                [-rho_1*rho_2*mass, MMOI[1], -rho_2*rho_3*mass],
                [-rho_1*rho_3*mass, -rho_2*rho_3*mass, MMOI[2]]]
        
<<<<<<< HEAD

=======
        
        omega_bar = np.matmul(omega_mat[minute], bees[minute])
        
        inter = np.matmul(J_mat, bees[minute])
        J_bar = np.matmul(np.transpose(bees[minute]), inter)
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
        
        # =============================================================================
        # The resulting equation in which the external moments and control moments need
        # to be added
        # =============================================================================
        
<<<<<<< HEAD
        omega_dot[0] = ((tau_build[minute][0]) + ((MMOI[1]-MMOI[2]) * omega_mat[minute][1] * omega_mat[minute][2]))/ MMOI[0]
        omega_dot[1] = ((tau_build[minute][1]) + ((MMOI[2]-MMOI[0]) * omega_mat[minute][2] * omega_mat[minute][0]))/ MMOI[1]
        omega_dot[2] = ((tau_build[minute][2]) + ((MMOI[0]-MMOI[1]) * omega_mat[minute][0] * omega_mat[minute][1]))/ MMOI[2]
=======
        omega_dot = omega_bar**4 - J_bar**-1 * (tau_orb[minute][0:3])
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
        
        attitude = [attitude[0] + (omega_mat[minute][0])*60,
                        attitude[1] + (omega_mat[minute][1])*60, attitude[2] + omega_mat[minute][2]*60]
        
        omega_mat.append([omega_mat[minute][0] + (omega_dot[0]), 
                          omega_mat[minute][1]+ (omega_dot[1]), omega_mat[minute][2] + (omega_dot[2])])
        
<<<<<<< HEAD
    print(attitude)
    return omega_mat[-1], minute



EOE_omega = calculate_angular_velocity()

print("Time to detumble after eclipse: ", calculate_detumble_time(15, EOE_omega, 7))
=======

    return omega_mat[-1], minute


for j in range(1,30):
    k = j/10
    EOE_omega, count = calculate_angular_velocity()
    print("Time to detumble after eclipse: ", calculate_detumble_time(k, EOE_omega, count))
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
