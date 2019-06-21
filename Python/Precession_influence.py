# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:08:45 2019

@author: tom-h_000
"""

import numpy as np
from Magnetorquer_calc import tau_orb_crop, get_Js, angle_magn
from ADCS_Mag import magn_P_data
import os
import matplotlib.pyplot as plt
from ADCS_combined import tau_orb, tau_build

# =============================================================================
# Variables
# =============================================================================

attitude = []

initial_attitude = [0, 0 , 0]

attitude.append(initial_attitude)
I_x = 0.013026
MMOI = [45.868,
       45.868,
       83.298] #kg * m^2

mass = 250 #kg

count = 0

margin = 0.5 * (np.pi/180)

# =============================================================================
# From matlab the values for b for the orbit will have to be imported to get
# the most accurate results possible
# =============================================================================

file = open('bees.txt', "r")
lines = file.readlines()                                                      #Read lines 
file.close()
b = []
for s in lines:                                 
    b.append(s.replace('i', '').replace('0-', '').replace('0+', '').replace(',', ' ').split())   

bees = np.asarray(b).astype(np.float)

# =============================================================================
# Calculate the angle of the sc with respect to the center of earth
# =============================================================================
def calculate_sc_angle(count, omegas):
    #using count-1 as count has not been written
    new_attitude = [attitude[count-1][0] + (omegas[0])*60,
                    attitude[count-1][1] + (omegas[1])*60, attitude[count-1][2] + omegas[2]*60]
    
    #append the new attitude to the lists with attitudes to be able to keep track
    for i in range(3):
        if not np.pi/2 < new_attitude[i] < np.pi/2:
            new_attitude[i] = 0
    attitude.append(new_attitude)
    
    phi_sc = attitude[count][2]
    omega_sc = attitude[count][1]
    theta_sc = attitude[count][0]
    
    return theta_sc, omega_sc, phi_sc
    
# =============================================================================
# Calculate the over all angle between the sc and the magnetic field
# =============================================================================
def calculate_angle(theta_B, phi_B, omega_B, theta_sc, phi_sc, omega_sc):
    theta = abs(theta_sc) - theta_B
    phi = abs(phi_sc) - phi_B
    omega = abs(omega_sc) - omega_B
    
    return theta, omega, phi
    

def calculate_detumble_time(initial_dipole, omega_matrix, minute = 0):
    count = 0;
    dipole = ([initial_dipole, initial_dipole, initial_dipole])
    # =============================================================================
    #Starting values for angular velocity    
    # =============================================================================
    omega_mat = omega_matrix
    
    
    
    while not (-margin < omega_mat[0] < margin) or not (-margin < omega_mat[1]< margin)  or not (-margin < omega_mat[2] < margin):
    
        for i in range(3):
            if omega_mat[i] > 0:
                dipole[i] = - abs(dipole[i])
            elif omega_mat[i] <= 0:
                dipole[i] = abs(dipole[i])
        
        orbit_torque = get_Js(dipole, minute)
    
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
        
        
        omega_bar = np.matmul(np.transpose(omega_mat), bees[minute])
        
        inter = np.matmul(J_mat, bees[minute])
        J_bar = np.matmul(np.transpose(bees[minute]), inter)
        
        # =============================================================================
        # The resulting equation in which the external moments and control moments need
        # to be added
        # =============================================================================
        
        omega_dot = omega_bar**4 - J_bar**-1 * (orbit_torque + (tau_build[minute]))
        
#        print(omega_dot)
#        print("angular velocities \n", omega_mat)
        
        omega_mat[0] = omega_mat[0] + ((omega_dot[0]))
        omega_mat[1] = omega_mat[1] + ((omega_dot[1]))
        omega_mat[2] = omega_mat[2] + ((omega_dot[2]))

        #update all used values, angles and other shizzle
        count = count + 1
        theta_B, omega_B, phi_B = angle_magn(tau_orb[count])
        theta_sc, omega_sc, phi_sc = calculate_sc_angle(count, omega_mat) 
        theta, omega, phi = calculate_angle(theta_B, phi_B, omega_B, theta_sc, phi_sc, omega_sc)
        
        dipole = ([(phi/(90*(np.pi/180)))*initial_dipole, (theta/(90*(np.pi/180)))*initial_dipole, (omega/(90*(np.pi/180)))*initial_dipole])

        
        if count == 430:
            return -2;

    return count


            
x = []
detumble_times1= []
detumble_times2= [] 
detumble_times3= []  
k = 1;
for k in range(k, 421, 1):
#    omega_mat1 = np.array([[1.5 * (np.pi/180)],
#          [0 * np.pi/180],
#          [1.75 * (np.pi/180)]]) #rad/s
    omega_mat2 = np.array([[0.75 * (np.pi/180)],
      [0.75 * np.pi/180],
      [1.75 * (np.pi/180)]]) #rad/s
#    omega_mat3 = np.array([[0 * (np.pi/180)],
#          [1.5* np.pi/180],
#          [1.75 * (np.pi/180)]]) #rad/s
    j = k/10
#    y1 = calculate_detumble_time(j, omega_mat1)
    y2 = calculate_detumble_time(j, omega_mat2)
#    y3 = calculate_detumble_time(j, omega_mat3)
#    detumble_times1.append(y1)
    detumble_times2.append(y2)
#    detumble_times3.append(y3)
    
    x.append(j)

plt.plot(x, detumble_times2, 'r')
#, x,  detumble_times2, 'b', x, detumble_times3, 'g')
plt.axvline(x=1)
plt.axvline(x=15)
plt.text(x[9], detumble_times2[9], '  ({}, {})'.format(x[9], detumble_times2[9]))
plt.text(x[149], detumble_times2[149]+ 8, '  ({}, {})'.format(x[149], detumble_times2[149]))
plt.title("Time to detumble for different values of linear dipole moment")
plt.xlabel("Linear dipole moment(Am^2)")
plt.ylabel("Time to detumble(min)")
plt.show()






    



