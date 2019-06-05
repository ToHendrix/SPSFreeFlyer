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


# =============================================================================
# Creating the Aerodynamic Database
# =============================================================================
def dataset(file):
    temp = []
    
    f = open(file,"r") # Open file for reading
    lines = f.readlines()
    f.close()
    for s in lines:
        temp.append(s.replace(',', '').split())
    
    if file == 'Aero_data.txt':
        del temp[0:89]
        del temp[-1]
        temp = np.asarray(temp).astype(np.float)
#        temp = np.vstack((np.array(['N_0 [cm^-3]', 'FFlux [cm^-2 s^-1]', 
#                                    'BFlux [cm^-2 s^-1]', 'FFluence [cm^-2]', 
#                                    'BFluence [cm^-2]', 'TMD [g cm^-3]', 
#                                    'Angle [deg]']), temp))
    if file == 'Coordinate_data.txt':
        del temp[0:73]
        del temp[-1]
        temp = np.asarray(temp).astype(np.float)
#        temp = np.vstack((np.array(['MJD [days]', 'Altitude [km]', 'Latitude [deg]', 
#                                    'Longitude [deg]', 'Local Time [hrs]', 
#                                    'Alpha [deg]']), temp))
    if file == 'Vector_data.txt':
        del temp[0:73]
        del temp[-1]
#       -Attitude contains the 9 values for rotation matrix w.r.t. reference 
#        frame (1_rho, 1_theta, 1_phi) (https://orfeo.kbr.be/bitstream/handle/internal/6594/Heynderickx%281996j%29.pdf?sequence=1&isAllowed=y)
#       -Axes order: , X (In orbital plane, away from earth), 
#        Y (South, perpendicular to orbital plane), Z (parallel to V)
        temp = np.asarray(temp).astype(np.float)
#        temp = np.vstack((np.array(['Attitude_0 [-]', 'Attitude_1 [-]', 
#                                    'Attitude_2 [-]', 'Attitude_3 [-]', 
#                                    'Attitude_4 [-]', 'Attitude_5 [-]', 
#                                    'Attitude_6 [-]', 'Attitude_7 [-]', 
#                                    'Attitude_8 [-]', 'V_x [km s^-1]', 
#                                    'V_y [km s^-1]', 'V_z [km s^-1]', 
#                                    'Sun_norm_x [-]', 'Sun_norm_y [-]', 
#                                    'Sun_norm_z [-]']), temp))
    if file == 'Magnetic_field_data.txt':
        del temp[0:88]
        del temp[-1]
#       -B_SPH_RHO is the projection of the magnetic field on the radial axis, 
#        i.e. the vertical component of the magnetic field;
#   	-B_SPH_THETA is the southward component of the magnetic field after 
#        projection in the horizontal plane;
#       -B_SPH_PHI is the eastward component of the magnetic field after 
#        projection in the horizontal plane
        temp = np.asarray(temp).astype(np.float)
#        temp = np.vstack((np.array(['B_LOC [Gauss]','L_LOC [?]','B_ALP [Gauss]'
#                                    ,'L_ALP [?]', 'I_LOC [?]','I_ALP [Re]',
#                                    'B_SPH_RHO [Gauss]', 'B_SPH_THETA [Gauss]',
#                                    'B_SPH_PHI [Gauss]','N_SPH_RHO [-?]',
#                                    'N_SPH_THETA [-?]','N_SPH_PHI [-?]']), 
#                                    temp))
    
    return temp

aero_data = dataset('Aero_data.txt')
coor_data = dataset('Coordinate_data.txt')
vect_data = dataset('Vector_data.txt')
magn_data = dataset('Magnetic_field_data.txt')