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
    if file == 'Coordinate_data.txt':
        del temp[0:73]
        del temp[-1]
    if file == 'Vector_data.txt':
        del temp[0:73]
        del temp[-1]
    if file == 'Magnetic_field_data.txt':
        del temp[0:88]
        del temp[-1]
    

    return temp

aero = dataset('Aero_data.txt')
coor = dataset('Coordinate_data.txt')
vect = dataset('Vector_data.txt')
magn = dataset('Magnetic_field_data.txt')