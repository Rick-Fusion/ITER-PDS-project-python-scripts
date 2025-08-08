import imaspy
import imaspy.util
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# This python script is dedicated to plotting contours of the poloidal magnetic flux (psi) in a 2D plane from the equilibrium IDS at a certain time slice




matplotlib.use("TkAgg")

entry = imaspy.DBEntry(imaspy.ids_defs.HDF5_BACKEND, "ITER", 666666, 0, "vanschr", data_version="4")	#  define input data entry

entry.open()            # open input data entry

time_requested = 225                  # time requested

#eq = entry.get_slice("equilibrium", time_requested, imaspy.ids_defs.CLOSEST_INTERP)
eq = entry.get("equilibrium")

time = eq.time[0]

print("time at which profile is plotted:", time)



# Plasma boundary points

plasma_boundary_r = eq.time_slice[0].boundary.outline.r

plasma_boundary_z = eq.time_slice[0].boundary.outline.z


# Poloidal magnetic flux (psi) data

#r = eq.time_slice[0].profiles_2d[0].r

#z = eq.time_slice[0].profiles_2d[0].z

#psi  = eq.time_slice[0].profiles_2d[0].psi


r = eq.time_slice[0].ggd[0].r[0].values

z = eq.time_slice[0].ggd[0].z[0].values

psi  = eq.time_slice[0].ggd[0].psi[0].values



# Select which values of psi you want to plot a contour line

psi_min = psi.min()
psi_max = psi.max()

print("Minimum of psi =", psi_min, "and Maximum of psi =", psi_max)

psi_value_plasma_boundary = eq.time_slice[0].boundary.psi 		# Define value of psi at plasma boundary / separatrix

contour_values_psi = np.arange(psi_min, psi_max, 10)

print("Contour psi values are", contour_values_psi)

# Plot of Equilibrium data

plt.figure(figsize=(8, 8))

# Plot of the seperatrix

plt.plot(plasma_boundary_r, plasma_boundary_z, '-', color="red")


# The following would not only plot the seperatrix
#plt.tricontour(r.flatten(), z.flatten(), psi.flatten(), levels = [psi_value_plasma_boundary] , linewidths=2, colors=["red"])

# Plot of the magnetic flux contours
plt.tricontour(r.flatten(), z.flatten(), psi.flatten(), levels = contour_values_psi, linewidths=1, cmap="viridis")
plt.colorbar(label="Poloidal magnetic flux (Wb)")

# Additional scatter plot of psi

plt.scatter(r, z, c=psi, cmap="viridis", alpha=0.1)


# Labels and title
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.xlim(1, 10)
plt.ylim(-6, 6)

plt.show()

