import imaspy
import imaspy.util
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# This python script is dedicated to plotting plasma boundary in a 2D plane from the equilibrium IDS at a certain time slice




matplotlib.use("TkAgg")

entry = imaspy.DBEntry(imaspy.ids_defs.HDF5_BACKEND, "ITER", 666666, 1, "vanschr", data_version="4")	#  define input data entry

entry.open()            # open input data entry

time_requested = 225                  # time requested

eq = entry.get_slice("equilibrium", time_requested, imaspy.ids_defs.CLOSEST_INTERP)

time = eq.time[0]

print("time at which profile is plotted:", time)

plasma_boundary_r = eq.time_slice[0].boundary.outline.r

plasma_boundary_z = eq.time_slice[0].boundary.outline.z


plt.scatter(plasma_boundary_r, plasma_boundary_z, marker='o', s = 3)
plt.title("Plasma boundary outline in 2D poloidal plane at t={time} seconds")
plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.xlim(0,10)
plt.ylim(-6,6)

plt.show()
