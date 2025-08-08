import imas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.lines import Line2D
matplotlib.use("TkAgg")

""" The aim of this plotting script is to visualise PF coil currents over time for multiple iterations resulting from the NICE inverse simulation in the NICE inverse + TORAX coupling scheme """




#open the Data Entry
entry_1 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_1/post_inverse","r")
entry_2 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_2/post_inverse","r")
entry_3 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_3/post_inverse","r")


# Load data from all entries
pf_active_1 = entry_1.get("pf_active") 
time_1 = pf_active_1.time

pf_active_2 = entry_2.get("pf_active") 
time_2 = pf_active_2.time

pf_active_3 = entry_3.get("pf_active") 
time_3 = pf_active_3.time



# For each iteration, create 2D array: each coil vs time 
Ic_1 = np.zeros((12, len(time_1)))
time_2d_1 = np.tile(time_1, (12, 1))

Ic_2 = np.zeros((12, len(time_2)))
time_2d_2 = np.tile(time_2, (12, 1))

Ic_3 = np.zeros((12, len(time_3)))
time_2d_3 = np.tile(time_3, (12, 1))



# For each iteration, set coil current data to 2D array
for i in range(12):
    Ic_1[i] = pf_active_1.coil[i].current.data.value
    Ic_2[i] = pf_active_2.coil[i].current.data.value
    Ic_3[i] = pf_active_3.coil[i].current.data.value



# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(12, 10))

ax1.set_title(rf"PF coil currents ($I_p$) against time for all iterations")
colors = plt.get_cmap('tab20') 

coil_legend = []
iter_legend = []
linestyles = ['-', '--', ':']
iteration_labels = ['it=1', 'it=2', 'it=3']
coil_names = ['CS3U', 'CS2U', 'CS1U', 'CS1L', 'CS2L', 'CS3L', 'PF1', 'PF2', 'PF3', 'PF4', 'PF5', 'PF6']     # short coil names


# Plotting all the PF coil current data for each iteration
for i in range(12):
    coil_name = coil_names[i]

    color = colors(i)

    ax1.plot(time_2d_1[i,:], Ic_1[i,:]/1e3, linestyle=linestyles[0], color=color)
    ax1.plot(time_2d_2[i,:], Ic_2[i,:]/1e3, linestyle=linestyles[1], color=color)
    ax1.plot(time_2d_3[i,:], Ic_2[i,:]/1e3, linestyle=linestyles[2], color=color)

    # For coil legend: only add once per coil
    coil_legend.append(Line2D([0], [0], color=color, lw=2, label=coil_name))


# For iteration legend (line style only)
for ls, label in zip(linestyles, iteration_labels):
    iter_legend.append(Line2D([0], [0], color='black', linestyle=ls, lw=2, label=label))

# Add legends
legend1 = ax1.legend(handles=coil_legend, title='Coils', bbox_to_anchor=(1, 1), loc='upper left')
legend2 = ax1.legend(handles=iter_legend, title='Iterations', bbox_to_anchor=(1, 0.65), loc='upper left')

# Formatting
#ax1.set_aspect('equal')
ax1.set_xlabel(r"time (s)")
ax1.set_ylabel(r"$I_c$ (kA)")
#ax1.set_xlim(0, 1)
#ax1.set_ylim(-25, -20)
ax1.relim()
ax1.autoscale_view(True, True, True)

# Add both legends to the same axes
ax1.add_artist(legend1)

plt.tight_layout()
plt.show()


entry_1.close()
entry_2.close()
entry_3.close()
