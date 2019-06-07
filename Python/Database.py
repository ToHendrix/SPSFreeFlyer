# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 13:33:13 2019

@author: J.J.M. Frijns
"""

# =============================================================================
# This program contains the time step data used for the models implemented
# in the design of the ADCS subsystem
# =============================================================================

# =============================================================================
# Importing used packages, systems and files
# =============================================================================
import numpy as np

# =============================================================================
# Defining arrays and variables
# =============================================================================
f       = 1./298.2                                                             #Flattening of the earth [-]
a_Earth = 6378160.                                                             #Semi-major axis Earth [m]
b_Earth = a_Earth*(1. - f)                                                     #Semi-minor axis Earth [m] 

# =============================================================================
# Creating the databases
# =============================================================================
def dataset(file):                                                             #Definition that reads all data files to arrays
    temp = []                                                                  #Temporary to append an array line to      
    f = open(file,"r")                                                         #Open file for reading
    lines = f.readlines()                                                      #Read lines 
    f.close()                                                                  #Close file         
    for s in lines:                                 
        temp.append(s.replace(',', '').split())                                #Append a line to temp of which commas have been removed and has been split with the help of spacebars     
    
    if file == 'Aero_data.txt':
        del temp[0:89]                                                         #Deleting all non-numbers in front 
        del temp[-1]                                                           #Deleting all non-numbers in back
        temp = np.asarray(temp).astype(np.float)                               #Make an array of the list and floats of all values inside
#        temp = np.vstack((np.array(['N_0 [cm^-3]', 'FFlux [cm^-2 s^-1]',      #Legend 
#                                    'BFlux [cm^-2 s^-1]', 'FFluence [cm^-2]', 
#                                    'BFluence [cm^-2]', 'TMD [g cm^-3]', 
#                                    'Angle [deg]']), temp))
    if file == 'Coordinate_data.txt':
        del temp[0:73]                                                         #Deleting all non-numbers in front              
        del temp[-1]                                                           #Deleting all non-numbers in back             
        temp = np.asarray(temp).astype(np.float)                               #Make an array of the list and floats of all values inside 
        temp[:,2] = np.deg2rad(temp[:,2])                                      #Change 3rd column from degrees to radians     
        temp[:,3] = np.deg2rad(temp[:,3])                                      #Change 4th column from degrees to radians      
#        temp = np.vstack((np.array(['MJD [days]', 'Altitude [km]', 'Latitude  #Legend 
#                                    [rad]', 'Longitude [rad]', 'Local Time 
#                                    [hrs]', 'Alpha [deg]']), temp))

    if file == 'Vector_data.txt':
        del temp[0:73]                                                         #Deleting all non-numbers in front    
        del temp[-1]                                                           #Deleting all non-numbers in back       
#       -Attitude contains the 9 values for rotation matrix w.r.t. reference 
#        frame (1_rho, 1_theta, 1_phi) (https://orfeo.kbr.be/bitstream/handle/internal/6594/Heynderickx%281996j%29.pdf?sequence=1&isAllowed=y)
#       -Axes order: , X (In orbital plane, away from earth), 
#        Y (South, perpendicular to orbital plane), Z (parallel to V)
        temp = np.asarray(temp).astype(np.float)                               #Make an array of the list and floats of all values inside      
#        temp = np.vstack((np.array(['Attitude_0 [-]', 'Attitude_1 [-]',       #Legend 
#                                    'Attitude_2 [-]', 'Attitude_3 [-]', 
#                                    'Attitude_4 [-]', 'Attitude_5 [-]', 
#                                    'Attitude_6 [-]', 'Attitude_7 [-]', 
#                                    'Attitude_8 [-]', 'V_x [km s^-1]', 
#                                    'V_y [km s^-1]', 'V_z [km s^-1]', 
#                                    'Sun_norm_x [-]', 'Sun_norm_y [-]', 
#                                    'Sun_norm_z [-]']), temp))
    if file == 'Magnetic_field_data.txt':
        del temp[0:88]                                                         #Deleting all non-numbers in front         
        del temp[-1]                                                           #Deleting all non-numbers in back              
#       -B_SPH_RHO is the projection of the magnetic field on the radial axis, 
#        i.e. the vertical component of the magnetic field;
#   	-B_SPH_THETA is the southward component of the magnetic field after 
#        projection in the horizontal plane;
#       -B_SPH_PHI is the eastward component of the magnetic field after 
#        projection in the horizontal plane
        temp = np.asarray(temp).astype(np.float)                               #Make an array of the list and floats of all values inside     
#        temp = np.vstack((np.array(['B_LOC [Gauss]','L_LOC [?]','B_ALP [Gauss]' #Legend
#                                    ,'L_ALP [?]', 'I_LOC [?]','I_ALP [Re]',
#                                    'B_SPH_RHO [Gauss]', 'B_SPH_THETA [Gauss]',
#                                    'B_SPH_PHI [Gauss]','N_SPH_RHO [-?]',
#                                    'N_SPH_THETA [-?]','N_SPH_PHI [-?]']), 
#                                    temp))   
    return temp                                                                #Returning the final dataset array 

aero_data = np.array(dataset('Aero_data.txt'))                                 #Calling up the creation of aero dataset 
coor_elli_data = np.array(dataset('Coordinate_data.txt'))                      #Calling up the creation of coordinate dataset  
vect_C_data = np.array(dataset('Vector_data.txt'))                             #Calling up the creation of vector dataset  
magn_data = np.array(dataset('Magnetic_field_data.txt'))                       #Calling up the creation of magnetic field dataset  

# =============================================================================
# Converting the coordinate database from ellipsoidal coordinates to Cartesian
# coordinates
# =============================================================================
def elli_to_C():                                                               #     
    temp = np.zeros((43201,5))
    for i in range(len(coor_elli_data)):
        N   = a_Earth*(1. - f*(2. - f) * \
                       (np.sin(coor_elli_data[i,2]))**2.)**-0.5                #The radius of curvature in the prime vertical [m]
        x_C	= (N + coor_elli_data[i,1]*10.**3.) * \
                       np.cos(coor_elli_data[i,2]) \
                       * np.cos(coor_elli_data[i,3])                           #x coordinate C frame [m]
        y_C	= (N + coor_elli_data[i,1]*10.**3.) * \
                       np.cos(coor_elli_data[i,2]) \
                       * np.sin(coor_elli_data[i,3])                           #y coordinate C frame [m]
        e   = (a_Earth**2. - b_Earth**2.)**0.5 / a_Earth                       #The first eccentricity [-]
        z_C	= (N*(1. - e**2.) + coor_elli_data[i,1]*10.**3.) * \
                       np.sin(coor_elli_data[i,2])                             #z coordinate C frame [m]
        temp[i] = [x_C, y_C, z_C, N, e]
    
    return temp

coor_C_data = elli_to_C()                                                      #x_C, y_C, z_C, N, e

# =============================================================================
# Converting the coordinate database from the C frame to the E frame
# =============================================================================
def C_to_E(C_frame, angles):                                                   #Definition to transform from the C- to the E-frame 
    temp = np.zeros((43201,3))                                                 #Creating initial array to append to 
    
    for i in range(len(C_frame[:,0])):
        T_EC = np.array([[-np.sin(angles[i,2])*np.cos(angles[i,3]), \
                          -np.sin(angles[i,3])*np.cos(angles[i,2]), \
                          np.cos(angles[i,2])],
                         [-np.sin(angles[i,3]), np.cos(angles[i,3]), 0],
                         [-np.cos(angles[i,2])*np.cos(angles[i,3]), \
                          -np.cos(angles[i,2])*np.sin(angles[i,3]), \
                          -np.sin(angles[i,2])]])
        
        temp[i] = np.dot(T_EC, C_frame[i,:3])                                  #Transforming every data point in the line used 
        
    return temp

coor_E_data = C_to_E(coor_C_data, coor_elli_data)                              #Calling up the transformation from coordinate C- to E-frame 

# =============================================================================
# Converting the coordinate database from the E frame to the P frame
# =============================================================================
Sun_E_data = C_to_E(vect_C_data[:,-3:], coor_elli_data)                        # Converting the Sun vector from the E frame to the P frame 

def E_to_P(E_frame):                                                           #Definition to transform from the C- to the E-frame with the z-axis pointing at the Sun 
    temp = np.zeros((43201,3))
    
    for i in range(len(E_frame[:,0])):
        alpha = np.arctan2(Sun_E_data[i,1], Sun_E_data[i,2]) -\
        np.arctan2(coor_E_data[i,1], coor_E_data[i,2])  
        beta = np.arctan2(Sun_E_data[i,0], Sun_E_data[i,2]) - \
        np.arctan2(coor_E_data[i,0], coor_E_data[i,2])

        T_PE = np.array([[np.cos(beta), np.sin(beta)*np.sin(alpha), \
                          np.sin(beta)*np.cos(beta)],
                 [0., np.cos(alpha), -np.sin(alpha)],
                 [-np.sin(beta), np.cos(beta)*np.sin(alpha), \
                  np.cos(beta)*np.cos(alpha)]])
        
        temp[i] = np.dot(T_PE, E_frame[i,:3])
        
    return temp

coor_P_data = E_to_P(coor_E_data)





































