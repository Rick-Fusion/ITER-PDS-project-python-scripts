import imaspy
import imaspy.util
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# This python script is dedicated to plotting poloidal magnetic flux (psi) in a 2D plane from the equilibrium IDS at a certain time slice




matplotlib.use("TkAgg")

entry = imaspy.DBEntry(imaspy.ids_defs.HDF5_BACKEND, "ITER", 666666, 0, "vanschr", data_version="4")	#  define input data entry

entry.open()            # open input data entry

time_requested = 225                  # time requested

#eq = entry.get_slice("equilibrium", time_requested, imaspy.ids_defs.CLOSEST_INTERP)
eq = entry.get("equilibrium")

time = eq.time[0]

print("time at which profile is plotted:", time)


plasma_boundary_r = eq.time_slice[0].boundary.outline.r

plasma_boundary_z = eq.time_slice[0].boundary.outline.z


#plt.scatter(plasma_boundary_r, plasma_boundary_z, marker='o', s = 3)
#plt.title("Plasma boundary outline in 2D poloidal plane at t={time} seconds")
#plt.xlabel("R (m)")
#plt.ylabel("Z (m)")
#plt.xlim(0,10)
#plt.ylim(-6,6)



# When 2d profile data is in profiles_2d
#r = eq.time_slice[0].profiles_2d[0].r

#z = eq.time_slice[0].profiles_2d[0].z

#psi  = eq.time_slice[0].profiles_2d[0].psi



# When 2d profile data is in ggd (for equilibria that are computed with direct mode)

r = eq.time_slice[0].ggd[0].r[0].values

z = eq.time_slice[0].ggd[0].z[0].values

psi  = eq.time_slice[0].ggd[0].psi[0].values








plt.figure(figsize=(8, 8))
scatter = plt.scatter(r, z, c=psi, cmap='viridis', s=4)  # `s` sets the size of the points

plt.scatter(plasma_boundary_r, plasma_boundary_z,color='red', marker='o', s = 3)

# Add a colorbar to indicate the mapping of psi_normalized values to colors
plt.colorbar(scatter, label='Poloidal magnetic flux (Wb)')

# Labels and title
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.xlim(1, 10)
plt.ylim(-6, 6)

plt.show()

