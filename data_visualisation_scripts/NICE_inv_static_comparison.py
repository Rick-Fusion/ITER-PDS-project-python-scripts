import imas
import imas.util
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patches as patches



matplotlib.use("TkAgg")

# This python script is dedicated to visualising output data from NICE static inverse mode. 
# The input (desired) and output (resulting) equilibrium data is plotted in the 2D poloidal plane at the given time
# The original current data at the given time is shown (if available for comparison) and the computed current data at the given time is shown for each coil



# Loading IDSs from the desired location (input and output)

run_input = 'nice_inv_test_0'
run_output = 'nice_inv_test_3'

input_entry = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_coupling_test_directory/pds/run/pds/ymmsl_files/output/sink_nice_inv_test_0","r")
output_entry = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_coupling_test_directory/pds/run/pds/ymmsl_files/output/sink_nice_inv_test_3","r")



# Get the input data only at the requested time slice for which static inverse mode is requested / computed
eq_full = input_entry.get("equilibrium", lazy=True)
time_array = eq_full.time
    

time_inverse_mode = time_array[75]                       # time requested inverse mode

eq_input = input_entry.get_slice("equilibrium", time_inverse_mode, imas.ids_defs.CLOSEST_INTERP)
pf_active_original = input_entry.get_slice("pf_active", time_inverse_mode, imas.ids_defs.CLOSEST_INTERP)

# Get the output data computed from inverse mode

eq_output = output_entry.get_slice("equilibrium", time_inverse_mode, imas.ids_defs.CLOSEST_INTERP)
pf_active_output = output_entry.get_slice("pf_active", time_inverse_mode, imas.ids_defs.CLOSEST_INTERP)


###
# Assign required input data (equilibrium data and coil currents)

number_of_coils = pf_active_original.coil.size
time = eq_input.time[0]                         # Actual time at which equilibrium is taken from data


# Plasma boundary points - to plot the seperatrix in the equilibrium contour plot

plasma_boundary_r = eq_input.time_slice[0].boundary.outline.r

plasma_boundary_z = eq_input.time_slice[0].boundary.outline.z


# Poloidal magnetic flux (psi) data


# In run 1 (DINA data) the poloidal magnetic flux data is stored under profiles_2D, this is the correct way of extracting this data
'''
r = eq_input.time_slice[0].profiles_2d[0].r

z = eq_input.time_slice[0].profiles_2d[0].z

psi  = eq_input.time_slice[0].profiles_2d[0].psi
'''
# In equilibrium IDSs computed with NICE inverse the poloidal magnetic flux data is stored under ggd, this is the correct way of extracting this data


r = eq_input.time_slice[0].ggd[0].r[0].values

z = eq_input.time_slice[0].ggd[0].z[0].values

psi = eq_input.time_slice[0].ggd[0].psi[0].values



# 'original' coil current data (put in np array) - (only applicable if this data exists)

coil_currents_original = np.zeros(number_of_coils)

for i in range(number_of_coils):
    coil_currents_original[i] = pf_active_original.coil[i].current.data[0]




###


# Assign required output data (equilibrium and coil currents)

# inverse equilibrium data (flux map and plasma contour)

# In run 1 (DINA data) or data coming from this, the poloidal magnetic flux data is stored under profiles_2D, this is the correct way of extracting this data


'''
r_output = eq_output.time_slice[0].profiles_2d[0].r

z_output = eq_output.time_slice[0].profiles_2d[0].z

psi_output  = eq_output.time_slice[0].profiles_2d[0].psi

'''

# In equilibrium IDSs computed with NICE inverse the poloidal magnetic flux data is stored under ggd, this is the correct way of extracting this data


r_output = eq_output.time_slice[0].ggd[0].r[0].values

z_output = eq_output.time_slice[0].ggd[0].z[0].values

psi_output  = eq_output.time_slice[0].ggd[0].psi[0].values





plasma_boundary_r_output = eq_output.time_slice[0].boundary.outline.r

plasma_boundary_z_output = eq_output.time_slice[0].boundary.outline.z


# inverse coil current data

coil_currents_output = np.zeros(number_of_coils)

for i in range(number_of_coils):
    coil_currents_output[i] = pf_active_output.coil[i].current.data[0]


# Select which values of psi you want to plot a contour line (want to find the minimum and maximum for both equilibrium to decide ultimate max and min)

psi_min_input = psi.min()
psi_max_input = psi.max()

psi_min_output = psi_output.min()
psi_max_output = psi_output.max()

psi_min = min(psi_min_input, psi_min_output)
psi_max = max(psi_max_input, psi_max_output)


contour_values_psi = np.arange(psi_min, psi_max, 10)



# Now make a plot of the data (plot input equilibrium (when exist) , output equilibrium (by inverse), comparison coil currents input (if exists) and output in one plot)

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(10,10))


# Start with the equilibrium data plot


# Function to add coil rectangles to a plot (since want to do it for both input as output equilibrium plot)
def add_coil_rectangles(ax, pf_active_original, number_of_coils):
    for i in range(number_of_coils):
        coil_middle_r = pf_active_original.coil[i].element[0].geometry.rectangle.r
        coil_middle_z = pf_active_original.coil[i].element[0].geometry.rectangle.z
        coil_width = pf_active_original.coil[i].element[0].geometry.rectangle.width
        coil_height = pf_active_original.coil[i].element[0].geometry.rectangle.height
        
        # Calculate bottom-left corner from center
        bottom_left_R = coil_middle_r - coil_width / 2
        bottom_left_Z = coil_middle_z - coil_height / 2
        
        # Create and add the rectangle
        rect = patches.Rectangle(
            (bottom_left_R, bottom_left_Z),  # Bottom-left corner
            coil_width,                       # Width
            coil_height,                      # Height
            edgecolor="black",                # Outline color
            facecolor="none",                 # Transparent fill
            linewidth=2                       # Outline thickness
        )
        ax.add_patch(rect)


add_coil_rectangles(ax1, pf_active_original, number_of_coils)
add_coil_rectangles(ax2, pf_active_original, number_of_coils)


# The last two coils, Upper and Lower vertical Stabilization (VSL and VSU)coils, do not have rectangular shape, but the outline must be taken

VSU_outline_r = pf_active_original.coil[12].element[0].geometry.outline.r
VSU_outline_z = pf_active_original.coil[12].element[0].geometry.outline.z

VSL_outline_r = pf_active_original.coil[13].element[0].geometry.outline.r
VSL_outline_z = pf_active_original.coil[13].element[0].geometry.outline.z


# Plotting VS coils for input and output equilibrium

ax1.plot(VSU_outline_r, VSU_outline_z, 'o', markersize=0.6, color="black")
ax1.plot(VSL_outline_r, VSL_outline_z, 'o', markersize=0.6, color="black")

ax2.plot(VSU_outline_r, VSU_outline_z, 'o', markersize=0.6, color="black")
ax2.plot(VSL_outline_r, VSL_outline_z, 'o', markersize=0.6, color="black")



# Plot of the seperatrix for input and output equilibrium

ax1.plot(plasma_boundary_r, plasma_boundary_z, '-', color="red")

ax2.plot(plasma_boundary_r_output, plasma_boundary_z_output, '-', color="red", alpha=0.5, label='Resulting') 
ax2.plot(plasma_boundary_r, plasma_boundary_z, '--', color="red", label='Desired')               # show also desired plasma boundary in the inverse result plot


# Plot of the magnetic flux contours for input and output equilibrium

ax1.tricontour(r.flatten(), z.flatten(), psi.flatten(), levels = contour_values_psi, linewidths=1, cmap="viridis_r")

ax2.tricontour(r_output.flatten(), z_output.flatten(), psi_output.flatten(), levels = contour_values_psi, linewidths=1, cmap="viridis_r")



# Define normalization for consistent color mapping
norm = mcolors.Normalize(vmin=psi_min, vmax=psi_max)

# Custom color scale
#norm = mcolors.Normalize(vmin=-100, vmax=100)

# Additional scatter plot of psi

ax1.scatter(r, z, c=psi, cmap="viridis_r", alpha=0.05, norm=norm)

ax2.scatter(r_output, z_output, c=psi_output, cmap="viridis_r", alpha=0.05, norm=norm)

# Create a hidden color mapping object for the color bar
sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
sm.set_array([])  # Required for the color bar


# Create the color bar on a separate axis, adjust its position using `cax`
cbar_ax = fig.add_axes([0.34, 0.15, 0.015, 0.7])  # [left, bottom, width, height]
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label("Poloidal magnetic flux (Wb)")



# Formatting the plot

ax1.set_aspect('equal')  # Equal aspect ratio
ax1.set_xlabel("R (m)")
ax1.set_ylabel("Z (m)")
#ax1.set_title(f"Desired equilibrium poloidal flux contours at time: {time:.2f} s + pf_active Coil Cross-Sections")
ax1.set_title(f"equilibrium poloidal flux from run: {run_input} at time: {time:.2f} s")
ax1.set_xlim(0, 13)  # Adjust limits based on your data
ax1.set_ylim(-10, 10)


ax2.set_aspect('equal')  # Equal aspect ratio
ax2.set_xlabel("R (m)")
ax2.set_ylabel("Z (m)")
ax2.set_title(f"Resulting equilibrium poloidal flux contours from Inverse at time: {time:.2f} s + pf_active Coil Cross-Sections")
ax2.set_title(f"equilibrium poloidal flux from run: {run_output} at time: {time:.2f} s")
ax2.set_xlim(0, 13)  # Adjust limits based on your data
ax2.set_ylim(-10, 10)
ax2.legend()


# Set up a list of coil names (manually here for simplicity)

coil_names = ['CS3U', 'CS2U', 'CS1U', 'CS1L', 'CS2L', 'CS3L',
              'PF1', 'PF2', 'PF3', 'PF4', 'PF5', 'PF6',
              'VSU', 'VSL']




# Making the coil current plot

ax3.scatter(coil_names, coil_currents_original/(1e3), color='red', s=20, label='Coil currents original')       # From original pf_active IDS
ax3.scatter(coil_names, coil_currents_output/(1e3), color='black', s=20, label='Coil currents output')       # From output pf_active IDS

ax3.set_xlabel('Coil Name')
ax3.set_ylabel('Current (kA)')
ax3.set_title(f'Coil Currents at time: {time:.2f}')
ax3.set_xticks(range(len(coil_names)))  # Set ticks at indices of coil_names
ax3.set_xticks([x - 0.5 for x in range(len(coil_names))], minor=True)  # Add gridlines between ticks
ax3.grid(which='minor', axis='x', linestyle='--', alpha=0.7)  # Gridlines between labels
ax3.grid(which='major', axis='y', linestyle='--', alpha=0.7)  # Normal gridlines for y-axis
ax3.legend()


# Adjust layout to prevent overlap
plt.subplots_adjust(right=0.9)

# Show the plots
plt.tight_layout()
plt.show()