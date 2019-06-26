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

<<<<<<< HEAD
simulation = 2
=======
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
attitude = []

initial_attitude = [0, 0 , 0]

attitude.append(initial_attitude)
<<<<<<< HEAD
=======
I_x = 0.013026
MMOI = [45.868,
       45.868,
       83.298] #kg * m^2
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7

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
    

<<<<<<< HEAD
def calculate_detumble_time(initial_dipole, omega_matrix, minute = 0, MMOI = [45.868, 45.868, 83.298], off = 1):
=======
def calculate_detumble_time(initial_dipole, omega_matrix, minute = 0):
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
    count = 0;
    dipole = ([initial_dipole, initial_dipole, initial_dipole])
    # =============================================================================
    #Starting values for angular velocity    
    # =============================================================================
    omega_mat = omega_matrix
<<<<<<< HEAD
    cycle = 0
    
    omega_dot = [0,0,0]                                                         #initial acceleration, assumed to be 0
=======
    
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
    
    
    while not (-margin < omega_mat[0] < margin) or not (-margin < omega_mat[1]< margin)  or not (-margin < omega_mat[2] < margin):
    
        for i in range(3):
            if omega_mat[i] > 0:
                dipole[i] = - abs(dipole[i])
            elif omega_mat[i] <= 0:
                dipole[i] = abs(dipole[i])
        
        orbit_torque = get_Js(dipole, minute)
<<<<<<< HEAD
        

        
        cycle = cycle + 1
=======
    
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
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
        
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
        omega_bar = np.matmul(np.transpose(omega_mat), bees[minute])
        
        inter = np.matmul(J_mat, bees[minute])
        J_bar = np.matmul(np.transpose(bees[minute]), inter)
<<<<<<< HEAD
        #new
        if (cycle < 11):
            omega_dot[0] = ((orbit_torque[0] *60 + tau_build[minute][0]) + ((MMOI[1]-MMOI[2]) * omega_mat[1] * omega_mat[2]))/ MMOI[0]
            omega_dot[1] = ((orbit_torque[1] *60 + tau_build[minute][1]) + ((MMOI[2]-MMOI[0]) * omega_mat[2] * omega_mat[0]))/ MMOI[1]
            omega_dot[2] = ((orbit_torque[2] *60 + tau_build[minute][2]) + ((MMOI[0]-MMOI[1]) * omega_mat[0] * omega_mat[1]))/ MMOI[2]
        
        elif (off > 1):
            for i in range(off-2):
                omega_dot[0] = ((tau_build[minute][0]) + ((MMOI[1]-MMOI[2]) * omega_mat[1] * omega_mat[2]))/ MMOI[0]
                omega_dot[1] = ((tau_build[minute][1]) + ((MMOI[2]-MMOI[0]) * omega_mat[2] * omega_mat[0]))/ MMOI[1]
                omega_dot[2] = ((tau_build[minute][2]) + ((MMOI[0]-MMOI[1]) * omega_mat[0] * omega_mat[1]))/ MMOI[2]
            
            
        elif (cycle == (10 + off)):
            omega_dot[0] = ((tau_build[minute][0]) + ((MMOI[1]-MMOI[2]) * omega_mat[1] * omega_mat[2]))/ MMOI[0]
            omega_dot[1] = ((tau_build[minute][1]) + ((MMOI[2]-MMOI[0]) * omega_mat[2] * omega_mat[0]))/ MMOI[1]
            omega_dot[2] = ((tau_build[minute][2]) + ((MMOI[0]-MMOI[1]) * omega_mat[0] * omega_mat[1]))/ MMOI[2]
            cycle = 0
=======
        
        # =============================================================================
        # The resulting equation in which the external moments and control moments need
        # to be added
        # =============================================================================
        
        omega_dot = omega_bar**4 - J_bar**-1 * (orbit_torque + (tau_build[minute]))
        
#        print(omega_dot)
#        print("angular velocities \n", omega_mat)
>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7
        
        omega_mat[0] = omega_mat[0] + ((omega_dot[0]))
        omega_mat[1] = omega_mat[1] + ((omega_dot[1]))
        omega_mat[2] = omega_mat[2] + ((omega_dot[2]))

        #update all used values, angles and other shizzle
        count = count + 1
<<<<<<< HEAD
        theta_B, omega_B, phi_B = angle_magn(tau_build[count])
        theta_sc, omega_sc, phi_sc = calculate_sc_angle(count, omega_mat) 
        theta, omega, phi = calculate_angle(theta_B, phi_B, omega_B, theta_sc, phi_sc, omega_sc)
        

#        dipole = ([(phi/(90*(np.pi/180)))*initial_dipole, (theta/(90*(np.pi/180)))*initial_dipole, (omega/(90*(np.pi/180)))*initial_dipole])

        dipole = ([initial_dipole, initial_dipole, initial_dipole])
        if count == 43200:
            return 100000;

    return count

# =============================================================================
# Simulation of different options for linear dipole
# =============================================================================
if simulation == 1:           
    x = []
    detumble_times1= []
    detumble_times2= [] 
    detumble_times3= []  
    k = 11;
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
        y2 = calculate_detumble_time(j, omega_mat2, 7)
    #    y3 = calculate_detumble_time(j, omega_mat3)
    #    detumble_times1.append(y1)
        detumble_times2.append(y2)
    #    detumble_times3.append(y3)
        
        x.append(j)
    
    plt.plot(x, detumble_times2, 'r')
    #, x,  detumble_times2, 'b', x, detumble_times3, 'g')
    plt.axvline(x=10)
    plt.axvline(x=15)
    plt.axvline(x=30)
    plt.text(x[89], detumble_times2[89]+ 100, '({}, {})'.format(x[89], detumble_times2[89]))
    plt.text(x[139], detumble_times2[139]+ 20, '({}, {})'.format(x[139], detumble_times2[139]))
    plt.text(x[289], detumble_times2[289]+ 20, '({}, {})'.format(x[289], detumble_times2[289]))
    plt.title("Time to detumble for different values of linear dipole moment")
    plt.xlabel("Linear dipole moment(Am^2)")
    plt.ylabel("Time to detumble(min)")
    plt.show()

# =============================================================================
# Simulation for different values of the magnatic field and other disturbance
# torques. So at a different time
# =============================================================================

#Make a list which will employ all values for the detumble time
if (simulation == 2):
    detumble_times_bar = []
    omega_mat = np.array([[0.75 * (np.pi/180)],
      [0.75 * np.pi/180],
      [1.75 * (np.pi/180)]]) #rad/s
    for start in range(1, 999):
        omega_mat = np.array([[0.75 * (np.pi/180)],
                               [0.75 * np.pi/180],
                               [1.75 * (np.pi/180)]]) #rad/s
        time = calculate_detumble_time(15, omega_mat, start)
        detumble_times_bar.append(time)
    
    for entry in range(len(detumble_times_bar) + 1):
        for i in range(1, len(detumble_times_bar)- entry):
            if detumble_times_bar[entry] > detumble_times_bar[entry+i] + i and not detumble_times_bar[entry+i] == -2:
                detumble_times_bar[entry] = detumble_times_bar[entry + i] + i
                
                
        
    plt.hist(detumble_times_bar, 75)
    plt.title("Variation of detumble time with changing magnetic field properties(n = 1000)")
    plt.ylabel("Amount of occurences (-)")
    plt.xlabel("Time to detumble(min)")
    plt.show()
# =============================================================================
# testing
# =============================================================================
if(simulation == 3):
    omega_mat = np.array([[0.75 * (np.pi/180)],
      [0.75 * np.pi/180],
      [1.75 * (np.pi/180)]]) #rad/s
    print(calculate_detumble_time(15, omega_mat, 0))

# =============================================================================
# Sensitivity analysis on MMOI and therefore configuration
# =============================================================================
if(simulation == 4):
    
    MMOI = [45.868, 45.868, 83.298]
    detumble_times_MMOIX = []
    detumble_times_MMOIY = []
    detumble_times_MMOIZ = []
    detumble_times_MMOIALL = []
    I_diff = 20
    x_var = []
    y_var = []
    z_var = []
    All_var = []
    for i in np.arange(MMOI[0] - I_diff, MMOI[0] + I_diff, 0.01):
        omega_mat = np.array([[0.75 * (np.pi/180)],
          [0.75 * np.pi/180],
          [1.75 * (np.pi/180)]]) #rad/s
        MMOI_update = i
        detumble_times_MMOIX.append(calculate_detumble_time(15, omega_mat, 7, (MMOI_update ,MMOI[1], MMOI[2])))
        x_var.append(i)
    for j in np.arange(MMOI[1] - I_diff, MMOI[1] + I_diff, 0.01):
        omega_mat = np.array([[0.75 * (np.pi/180)],
          [0.75 * np.pi/180],
          [1.75 * (np.pi/180)]]) #rad/s
        MMOI_updatedy = j
        detumble_times_MMOIY.append(calculate_detumble_time(15, omega_mat, 7, (MMOI[0],MMOI_updatedy, MMOI[2])))
        y_var.append(j)
    for k in np.arange(MMOI[2] - I_diff, MMOI[2] + I_diff, 0.01):
        omega_mat = np.array([[0.75 * (np.pi/180)],
          [0.75 * np.pi/180],
          [1.75 * (np.pi/180)]]) #rad/s
        MMOI_updatedz = k
        detumble_times_MMOIZ.append(calculate_detumble_time(15, omega_mat, 7, (MMOI[0],MMOI[1], MMOI_updatedz)))
        z_var.append(k)
    
        
    
    plt.plot(x_var, detumble_times_MMOIX, label = 'MMOI_x')
    plt.plot(y_var, detumble_times_MMOIY, label = 'MMOI_y')
    plt.plot(z_var, detumble_times_MMOIZ, label = 'MMOI_z')
    plt.title("Variation of detumble time with MMOI values")
    plt.ylabel("Detumble time (min)")
    plt.xlabel("MMOI(kg * m^2)")
    plt.legend(loc = 'upper left')
    plt.show()
    
if (simulation == 5):
        
    
    
    
=======
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


>>>>>>> eed9440db50205d29f7a782ccadc133f098239b7




    



