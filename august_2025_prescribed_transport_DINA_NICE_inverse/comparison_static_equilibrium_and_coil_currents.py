import imas
import imas.util
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patches as patches
import sys



matplotlib.use("TkAgg")


# This python script is dedicated to comparing an existing static magnetic equilibrium of a tokamak plasma and the corresponding coil currents at a requested time with a newly computed equilibrium + coil currents with an inverse simulation.
# This script can for example be used to compare DINA results with NICE inverse results. This is useful for PDS deliverable 2 for the case of prescribed transport. 
#  - contours of the poloidal magnetic flux (psi) in a 2D plane at requested time (for input data and output inverse data)
#  - poloidal magnetic flux map at requested time (for input data and output inverse data)
#  - PF and CS coil cross-sections
#  - Vacuum vessel outline
#  - limiter wall + divertor outline
#  - PF and CS coil currents at the requested time


# For the plotting titles specify where input and output are coming from (Most likey input from DINA and output from NICE inverse)
input_data_source = "DINA"
# input_data_source = "NICE inverse"
output_data_source = "NICE inverse"



# Specify if the input and output magnetic equilibrium data is stored under profiles_2d (DINA) or ggd (NICE inverse)
equilibrium_in_data_stored = 'ggd'      # either 'profiles_2d'  (DINA results) or  'ggd'  (NICE inverse results)

equilibrium_out_data_stored = 'ggd'      # either 'profiles_2d'  (DINA results) or  'ggd'  (NICE inverse results)



if equilibrium_in_data_stored not in ('profiles_2d', 'ggd'):
    print(f"equilibrium_in_data_stored = {equilibrium_in_data_stored} is not a valid option. Change it to a valid option for acessing equilibrim data.")
    sys.exit()

if equilibrium_out_data_stored not in ('profiles_2d', 'ggd'):
    print(f"equilibrium_out_data_stored = {equilibrium_out_data_stored} is not a valid option. Change it to a valid option for acessing equilibrim data.")
    sys.exit()
    

print(f"equilibrium_in_data_stored = {equilibrium_in_data_stored}")

print(f"equilibrium_out_data_stored = {equilibrium_out_data_stored}")



# Loading IDSs from the desired location (input and output)

entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/NICE_inverse_correct_boundary/flattop_125s","r")
entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/NICE_inverse_correct_boundary_correct_profiles/flattop_125s","r")


entry_machine_description = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/3","r") # TODO: take geometry information from a general source? Machine description database?


# Only if length of time array from output equilibrium data is greater than one, provide requested time value.
requested_time_value = 8


# Either the inverse data only contains a single time slice or specify the time when there are more available

eq_out = entry_out.get("equilibrium", lazy=True)
time_out_array = eq_out.time

if len(time_out_array) == 1:
    time_out_single_time_slice = time_out_array[0]

elif len(time_out_array) == 0:
    print("There is no time in equilibrium.time array")
    sys.exit()

else:
    time_out_single_time_slice = requested_time_value
    


# Load the input equilibrium and pf active at the single time slice

eq_in = entry_in.get_slice("equilibrium", time_out_single_time_slice, imas.ids_defs.CLOSEST_INTERP)
pf_active_in = entry_in.get_slice("pf_active", time_out_single_time_slice, imas.ids_defs.CLOSEST_INTERP)


# Load the output equilibrium (only necessary if there are more time slices, since already lazy loaded) and pf active at the single time slice

if len(time_out_array) > 1:
    eq_out = entry_out.get_slice("equilibrium", time_out_single_time_slice, imas.ids_defs.CLOSEST_INTERP)

pf_active_out = entry_out.get_slice("pf_active", time_out_single_time_slice, imas.ids_defs.CLOSEST_INTERP)


wall = entry_machine_description.get("wall") 


### Assign required input and output data (equilibrium data and coil currents) ###


number_of_coils_in = pf_active_in.coil.size
number_of_coils_out = pf_active_out.coil.size

print(f"Amount of coils in input data = {number_of_coils_in}")
print(f"Amount of coils in output data = {number_of_coils_out}")




# Plasma boundary points of input equilibrium data - to plot the seperatrix in the input equilibrium contour plot

plasma_boundary_r_in = eq_in.time_slice[0].boundary.outline.r

plasma_boundary_z_in = eq_in.time_slice[0].boundary.outline.z


# Plasma boundary points of input equilibrium data - to plot the seperatrix in the input equilibrium contour plot

plasma_boundary_r_out = eq_out.time_slice[0].boundary.outline.r

plasma_boundary_z_out = eq_out.time_slice[0].boundary.outline.z




# Poloidal magnetic flux (psi) data from input equilibrium


# When equilibrium data is under profiles_2d field: e.g. for DINA results

if equilibrium_in_data_stored == 'profiles_2d':
    r_in = eq_in.time_slice[0].profiles_2d[0].r
    z_in = eq_in.time_slice[0].profiles_2d[0].z
    psi_in = eq_in.time_slice[0].profiles_2d[0].psi

# when equilibrium data is under ggd field: e.g. for NICE inverse results
if equilibrium_in_data_stored == 'ggd':
    r_in = eq_in.time_slice[0].ggd[0].r[0].values
    z_in = eq_in.time_slice[0].ggd[0].z[0].values
    psi_in = eq_in.time_slice[0].ggd[0].psi[0].values



# Poloidal magnetic flux (psi) data from output equilibrium


# When equilibrium data is under profiles_2d field: e.g. for DINA results

if equilibrium_out_data_stored == 'profiles_2d':
    r_out = eq_out.time_slice[0].profiles_2d[0].r
    z_out = eq_out.time_slice[0].profiles_2d[0].z
    psi_out = eq_out.time_slice[0].profiles_2d[0].psi

# when equilibrium data is under ggd field: e.g. for NICE inverse results
if equilibrium_out_data_stored == 'ggd':
    r_out = eq_out.time_slice[0].ggd[0].r[0].values
    z_out = eq_out.time_slice[0].ggd[0].z[0].values
    psi_out = eq_out.time_slice[0].ggd[0].psi[0].values






# input and output coil current data (put in np array)

coil_currents_in = np.zeros(number_of_coils_in)
coil_currents_out = np.zeros(number_of_coils_out)


for i in range(number_of_coils_in):
    coil_currents_in[i] = pf_active_in.coil[i].current.data[0]

for i in range(number_of_coils_out):
    coil_currents_out[i] = pf_active_out.coil[i].current.data[0]





# Select which values of psi you want to plot a contour line (want to find the minimum and maximum for both equilibrium to decide ultimate max and min)


psi_min_global = min(psi_in.min(), psi_out.min())
psi_max_global = max(psi_in.max(), psi_out.max())


number_of_contours = 100                                    # Here specify the amount of psi iso contours that are visualised in the plot
contour_values_psi = np.linspace(psi_min_global, psi_max_global, number_of_contours) 


# Interesting value of the poloidal magnetic flux (psi) are on the magnetic axis and plasma boundary - these will be displayed below the plots for the input and output data

psi_axis_in = eq_in.time_slice[0].global_quantities.psi_axis.value
psi_boundary_in = eq_in.time_slice[0].global_quantities.psi_boundary.value

psi_axis_out = eq_out.time_slice[0].global_quantities.psi_axis.value
psi_boundary_out = eq_out.time_slice[0].global_quantities.psi_boundary.value






# Now make a plot of the data

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(12,10))


# equilibrium input and output data plot


# Function to add coil rectangles to a plot (since want to do it for both input as output equilibrium plot)
def add_coil_rectangles(ax, pf_active, number_of_coils):
    for i in range(number_of_coils):
        coil_middle_r = pf_active.coil[i].element[0].geometry.rectangle.r
        coil_middle_z = pf_active.coil[i].element[0].geometry.rectangle.z
        coil_width = pf_active.coil[i].element[0].geometry.rectangle.width
        coil_height = pf_active.coil[i].element[0].geometry.rectangle.height
        
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


add_coil_rectangles(ax1, pf_active_in, number_of_coils_in)
add_coil_rectangles(ax2, pf_active_out, number_of_coils_out)


# The last two coils, Upper and Lower vertical Stabilization (VSL and VSU)coils, do not have rectangular shape, but the outline must be taken

VSU_outline_r = pf_active_in.coil[12].element[0].geometry.outline.r
VSU_outline_z = pf_active_in.coil[12].element[0].geometry.outline.z

VSL_outline_r = pf_active_in.coil[13].element[0].geometry.outline.r
VSL_outline_z = pf_active_in.coil[13].element[0].geometry.outline.z


# Plotting VS coils for input and output equilibrium

ax1.plot(VSU_outline_r, VSU_outline_z, 'o', markersize=0.6, color="black")
ax1.plot(VSL_outline_r, VSL_outline_z, 'o', markersize=0.6, color="black")

ax2.plot(VSU_outline_r, VSU_outline_z, 'o', markersize=0.6, color="black")
ax2.plot(VSL_outline_r, VSL_outline_z, 'o', markersize=0.6, color="black")



# Plotting ITER inner and outer wall outline

r_wall_inner = wall.description_2d[0].vessel.unit[0].annular.centreline.r
z_wall_inner = wall.description_2d[0].vessel.unit[0].annular.centreline.z

r_wall_outer = wall.description_2d[0].vessel.unit[1].annular.centreline.r
z_wall_outer = wall.description_2d[0].vessel.unit[1].annular.centreline.z

ax1.plot(r_wall_inner, z_wall_inner, color="black", alpha=1)
ax1.plot(r_wall_outer, z_wall_outer, color="black", alpha=1)

ax2.plot(r_wall_inner, z_wall_inner, color="black", alpha=1)
ax2.plot(r_wall_outer, z_wall_outer, color="black", alpha=1)




# Plotting ITER divertor contour

r_limiter = wall.description_2d[0].limiter.unit[0].outline.r
z_limiter = wall.description_2d[0].limiter.unit[0].outline.z

r_divertor = wall.description_2d[0].limiter.unit[1].outline.r
z_divertor = wall.description_2d[0].limiter.unit[1].outline.z

ax1.plot(r_limiter, z_limiter, color="black", alpha=1)
ax1.plot(r_divertor, z_divertor, color="black", alpha=1)

ax2.plot(r_limiter, z_limiter, color="black", alpha=1)
ax2.plot(r_divertor, z_divertor, color="black", alpha=1)





# Plot of the seperatrix for input and output equilibrium

ax1.plot(plasma_boundary_r_in, plasma_boundary_z_in, '.', color="red", alpha=0.8, markersize=2)

ax2.plot(plasma_boundary_r_out, plasma_boundary_z_out, '-', color="red", alpha=0.8, label='Inverse boundary') 
ax2.plot(plasma_boundary_r_in, plasma_boundary_z_in, '.', color="blue", alpha=0.5, label='Desired boundary', markersize=1.5)            # show also desired plasma boundary in the inverse result plot


# Plot of the magnetic flux contours for input and output equilibrium

ax1.tricontour(r_in.flatten(), z_in.flatten(), psi_in.flatten(), levels = contour_values_psi, linewidths=1, cmap="viridis_r")

ax2.tricontour(r_out.flatten(), z_out.flatten(), psi_out.flatten(), levels = contour_values_psi, linewidths=1, cmap="viridis_r")



# Define normalization for consistent color mapping
norm = mcolors.Normalize(vmin=psi_min_global, vmax=psi_max_global)

# Custom color scale
#norm = mcolors.Normalize(vmin=-100, vmax=100)

# Additional scatter plot of psi

opaqueness_of_psi_scatter = 0

ax1.scatter(r_in, z_in, c=psi_in, cmap="viridis_r", alpha=opaqueness_of_psi_scatter, norm=norm)

ax2.scatter(r_out, z_out, c=psi_out, cmap="viridis_r", alpha=opaqueness_of_psi_scatter, norm=norm)

# Create a hidden color mapping object for the color bar
sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
sm.set_array([])  # Required for the color bar


# Create the color bar on a separate axis, adjust its position using `cax`
cbar_ax = fig.add_axes([0.315, 0.15, 0.015, 0.7])  # [left, bottom, width, height]
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label("Poloidal magnetic flux (Wb)")



# Formatting the plot

ax1.set_aspect('equal')  # Equal aspect ratio
ax1.set_xlabel("R (m)")
ax1.set_ylabel("Z (m)")
ax1.set_title(f"Input equilibrium poloidal magnetic flux from {input_data_source} at time: {time_out_single_time_slice:.4f} s", fontsize=12, fontweight="bold")
ax1.set_xlim(0, 13)  # Adjust limits based on your data
ax1.set_ylim(-10, 10)


# Display the poloidal magnetic flux (psi) at axis and boundary below the plot

fig.text(
    0.5, -0.07,
    rf"Input:  $ψ_{{axis}}$={psi_axis_in:.3e} Wb,  $ψ_{{boundary}}$={psi_boundary_in:.4e} Wb",
    ha="center", va="top", transform=ax1.transAxes, fontsize=12, fontweight="bold"
)


ax2.set_aspect('equal')  # Equal aspect ratio
ax2.set_xlabel("R (m)")
ax2.set_ylabel("Z (m)")
ax2.set_title(f"Output equilibrium poloidal magnetic flux from {output_data_source} at time: {time_out_single_time_slice:.4f} s", fontsize=12, fontweight="bold")
ax2.set_xlim(0, 13)  # Adjust limits based on your data
ax2.set_ylim(-10, 10)
ax2.legend()


# Display the poloidal magnetic flux (psi) at axis and boundary below the plot

fig.text(
    0.5, -0.07,
    rf"Output:  $ψ_{{axis}}$={psi_axis_out:.3e} Wb,  $ψ_{{boundary}}$={psi_boundary_out:.4e} Wb",
    ha="center", va="top", transform=ax2.transAxes, fontsize=12, fontweight="bold"
)



# Set up a list of coil names (manually here for simplicity)

coil_names = ['CS3U', 'CS2U', 'CS1U', 'CS1L', 'CS2L', 'CS3L',
              'PF1', 'PF2', 'PF3', 'PF4', 'PF5', 'PF6',
              'VSU', 'VSL']




# Making the coil current plot

ax3.scatter(coil_names, coil_currents_in/(1e3), color='red', s=20, label=f'Coil currents input (from {input_data_source})')       # From input pf_active IDS
ax3.scatter(coil_names, coil_currents_out/(1e3), color='black', s=20, label=f'Coil currents output (from {output_data_source})')       # From output pf_active IDS

ax3.set_xlabel('Coil Name')
ax3.set_ylabel('Current (kA)')
ax3.set_title(f'in- and ouput coil currents at time: {time_out_single_time_slice:.4f} s', fontsize=12, fontweight="bold")
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