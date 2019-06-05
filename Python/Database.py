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


# =============================================================================
# Defining arrays and variables
# =============================================================================
aero = []
coor = []
vect = []
magn = []

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
        temp.insert(0,['Density [cm^-3]', 'FFlux [cm^-2 s^-1]', 
                    'BFlux [cm^-2 s^-1]', 'FFluence [cm^-2]', 
                    'BFluence [cm^-2]', 'TMD [g cm^-3]', 'Angle [deg]'])
    if file == 'Coordinate_data.txt':
        del temp[0:73]
        del temp[-1]
        temp.insert(0,['MJD [days]', 'Altitude [km]', 'Latitude [deg]', 
                       'Longitude [deg]', 'Local Time [hrs]', 'Alpha [deg]'])
    if file == 'Vector_data.txt':
        del temp[0:73]
        del temp[-1]
        temp.insert(0,['Attitude_0 [?]', 'Attitude_1 [?]', 'Attitude_2 [?]', 
                       'Attitude_3 [?]', 'Attitude_4 [?]', 'Attitude_5 [?]', 
                       'Attitude_6 [?]', 'Attitude_7 [?]', 'Attitude_8 [?]', 
                       'V_1 [km s^-1]', 'V_2 [km s^-1]', 'V_3 [km s^-1]', 
                       'V_Sun_1 [km s^-1]', 'V_Sun_2 [km s^-1]', 
                       'V_Sun_3 [km s^-1]'])
    if file == 'Magnetic_field_data.txt':
        del temp[0:88]
        del temp[-1]
        temp.insert(0,['B_LOC [Gauss]','L_LOC [?]','B_ALP [Gauss]','L_ALP [?]',
                       'I_LOC [?]','I_ALP [Re]','B_SPH_RHO [Gauss]',
                       'B_SPH_THETA [Gauss]','B_SPH_PHI [Gauss]','N_SPH_1 [-?]',
                       'N_SPH_2 [-?]','N_SPH_3 [-?]'])
    

    return temp

aero = dataset('Aero_data.txt')
coor = dataset('Coordinate_data.txt')
vect = dataset('Vector_data.txt')
magn = dataset('Magnetic_field_data.txt')