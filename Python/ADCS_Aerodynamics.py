# Import necessary things
import scipy
from Database import *
#-----------------------------------------------------------------------------------------------------------------------
# ------------------------------------ Import Stuff from the Database --------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# velocity data from Database
Vx = 1000*vect_P_data[:,0] # X component of velocity in P-frame [m/s]
Vy = 1000*vect_P_data[:,1] # Y component of velocity in P-frame [m/s]
Vz = 1000*vect_P_data[:,2] # Z component of velocity in P-frame [m/s]
V = ((Vx)**2 + (Vy)**2 + (Vz)**2)**0.5 # Overall scalar velocity [m/s]
Vu = np.array([Vx, Vy, Vz])/V  # Unit velocity direction vector [m/s]

# Define Dimensions
r_plate = np.array([0, 0, .5 * h]) # This is the distance vector from cylinder coordinate system center to the plate center of pressure [m]
r_wall = np.array([Vu[0][:]*d/2, Vu[1][:]*d/2, np.zeros(43201)]) # Distance vector from center of cylinder coordinate system to outer center of pressure [m]
r_com = np.array([0.0015, 0.0015, 0*h/2]) # This is the vector CG location with respect to the geometric center [m]

# Density and Beta
beta = np.cos(-Vu[2])  # Angle defined in Oomen for P-frame F_wall decomposition [rad]
rho = 1000*aero_data[:,5] # Total Mass Density [kg/m3] from the Datasheet.

#-----------------------------------------------------------------------------------------------------------------------
# ------------------------------------ Define functions -----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def forces(C_D, rho, d, V, Vu, beta): # this computes and decomposes forces into P axes for later moment computation.
    F_plate = .5 * C_D * rho * (np.pi * (d / 2.)**2) * V**2 * (Vu) # F_plate is the force on the top of the cylinder (acting thru z line), as defined in Oomen [N]
    F_wall = -0.5 * C_D * rho * (2. * np.pi * d / 2. * h) * V**2 * (Vu * np.sin(beta)) #F_wall is force vector on the side of the cylinder (acting thru xy plane) as defined in Oomen [N]
    return F_plate, F_wall # outputs are the plate and wall forces, as defined in prev. lines.


def moments(r_plate, r_wall, F_plate, F_wall): # calculates the moments due to plate and wall forces in Nm, and sums to get the total moment due to aero M_ext_aero
    M_plate = [] # empty list for appending moment values due to difficulty of non-looped cross product
    M_wall = [] # same as prev.
    for i in range(0, 43201): # need to compute both moments for every point in the orbit (different velocity, attitude, density)
        M_plate.append(np.cross((r_plate-r_com), F_plate[:, i])) # append computed moment value. The formula is from Oomen again, and is just a cross product between the r vectors and the forces.
        M_wall.append(np.cross(((r_wall[:, i])-r_com), F_wall[:,i])) # note the r_wall-r_com: we need the distance between the CG and COP so that's what this gets.
    M_ext_aero = np.array(M_plate) + np.array(M_wall) # simply sum the wall and plate moments to get total M_ext_aero
    return M_ext_aero # outputs as promised. The shape for M_ext_aero is [[43201 x moment components], [43201 y moment components], [43201 z moment components]]


# ------------------------------------ Test (delet if annoying) -----------------------------------------------------------------------------
# Forces = forces(C_D, rho, d)
#print(Forces[1].shape)
# F_p = Forces[0]
# F_w = Forces[1]
# print(F_p[:,2].shape)
# print(F_w.shape)
# print(F_p[:,2]-F_p[:,1])
#print(F_p)

# Moments = moments(r_plate, r_wall, F_p, F_w)
# print(Moments[2][1])

#print(((vect_C_data[9][10]*1000)**2 + (vect_C_data[10][10]*1000)**2 + (vect_C_data[11][10]*1000)**2))
#
# print(len(aero_data[:,1]))
# print((aero_data[:,1]))
# print((aero_data[:,1])**2**0.5)
# print((r_plate.transpose()-r_com.transpose()).shape)
# print((r_plate-r_com))
# print((np.matrix(Forces[0][:,1])).shape)
# print((np.matrix(Forces[0][:,1])))
# print( np.cross((r_plate.transpose()-r_com.transpose()), Forces[0][:,1])[2] )
