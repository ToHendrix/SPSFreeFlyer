# Import necessary things
import scipy
import numpy as np
from Database import *

# ------------------------------------ Import Stuff from the Database --------------------------------------------------
# velocity data
Vx = vect_C_data[:,9] # km/s
Vy = vect_C_data[:,10] # km/s
Vz = vect_C_data[:,11] # km/s
V = ((1000*Vx)**2 + (1000*Vy)**2 + (1000*Vz)**2)**0.5 # m/s



# ------------------------------------ Define constants -----------------------------------------------------------------
h = 0.455  # [m], height of the SPS
d = 1.  # [m], diameter of SPS
C_D = 2.

# ------------------------------------ Define Vectors -------------------------------------------------------------------
Vu = np.array([Vx, Vy, Vz])/V  # [m/s], unit velocity vector
r_plate = np.array((0, 0, .5 * h))  # [m]
r_wall = np.array((0, 0, 9.11))  # [m] unsure about this, evaluate worst-case?

# ------------------------------------ assumed values -------------------------------------------------------------------
beta = 0.5  # [rad], this is wrong
rho = 0.0000000000001  # [kg/m3]


# ------------------------------------ Define functions -----------------------------------------------------------------
def forces(C_D, rho, d):
    F_plate = .5 * C_D * rho * (np.pi * d / 2.) * V**2 * (Vu)
    F_wall = -0.5 * C_D * (2. * np.pi * d / 2. * h) * V**2 * (Vu * np.sin(beta))
    return F_plate, F_wall


def moments(r_plate, r_wall, F_plate, F_wall):
    M_plate = np.cross(r_plate, F_plate)
    M_wall = np.cross(r_wall, F_wall)
    return M_plate, M_wall


# ------------------------------------ Test -----------------------------------------------------------------------------
Forces = forces(C_D, rho, d)
Moments = moments(r_plate, r_wall, Forces[0], Forces[1])
#print(((vect_C_data[9][10]*1000)**2 + (vect_C_data[10][10]*1000)**2 + (vect_C_data[11][10]*1000)**2))
#
# print(len(aero_data[:,1]))
# print((aero_data[:,1]))
# print((aero_data[:,1])**2**0.5)
print("a")