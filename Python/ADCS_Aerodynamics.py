# Import necessary things
import scipy
import numpy as np

# ------------------------------------ Define constants -----------------------------------------------------------------
h = 0.455  # [m], height of the SPS
d = 1.  # [m], diameter of SPS
C_D = 2.

# ------------------------------------ Define Vectors -------------------------------------------------------------------
Vu = np.array((1, 0, 0))  # [m/s], unit velocity vector. This is not yet determined, and needs a direction function
r_plate = np.array((0, 0, .5 * h))  # [m]
r_wall = np.array((0, 0, 9.11))  # [m] ("not", "sure", "yet")

# ------------------------------------ assumed values -------------------------------------------------------------------
V = 7000.  # [m/s]
beta = 0.5  # [rad], this is wrong
rho = 0.0000000000001  # [kg/m3]


# ------------------------------------ Define functions -----------------------------------------------------------------
def forces(C_D, rho, d):
    F_plate = .5 * C_D * rho * (np.pi * d / 2.) * V ** (2.) * (Vu)
    F_wall = -0.5 * C_D * (2. * np.pi * d / 2. * h) * V ** (2.) * (Vu * np.sin(beta))
    return F_plate, F_wall


def moments(r_plate, r_wall, F_plate, F_wall):
    M_plate = np.cross(r_plate, F_plate)
    M_wall = np.cross(r_wall, F_wall)
    return M_plate, M_wall


# ------------------------------------ Test -----------------------------------------------------------------------------
Forces = forces(C_D, rho, d)
Moments = moments(r_plate, r_wall, Forces[0], Forces[1])
print(Forces)
print(Moments)