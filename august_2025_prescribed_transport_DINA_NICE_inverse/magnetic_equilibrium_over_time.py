import imas
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from matplotlib.widgets import Slider
import sys


matplotlib.use("TkAgg")

# function for controlling slider with left and right keys

def on_key(event):
    """Handle key presses to adjust the slider value."""
    if event.key == 'right':
        # Increase slider value by 1 (move right)
        time_slider.set_val(min(time_slider.val + 1, time_slider.valmax))
    elif event.key == 'left':
        # Decrease slider value by 1 (move left)
        time_slider.set_val(max(time_slider.val - 1, time_slider.valmin))



# This python script is dedicated to plotting the magnetic equilibrium of a tokamak plasma, with possible time dependence. 
# This includes:
#  - contours of the poloidal magnetic flux (psi) in a 2D plane over time
#  - poloidal magnetic flux map over time
#  - PF and CS coil cross-sections
#  - Vacuum vessel outline
#  - limiter wall + divertor outline




# Specify if the equilibrium data is stored under profiles_2d (DINA) or ggd (NICE inverse)
equilibrium_data_stored = 'profiles_2d'      # either 'profiles_2d'  (DINA results) or  'ggd'  (NICE inverse results)

if equilibrium_data_stored not in ('profiles_2d', 'ggd'):
    print(f"equilibrium_data_stored = {equilibrium_data_stored} is not a valid option. Change it to a valid option for acessing equilibrim data.")
    sys.exit()
    

print(f"equilibrium_data_stored = {equilibrium_data_stored}")



# Open entries to read IDSs from
entry = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/DINA_data_3_105084_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_replaced_psi_boundary","r")
# entry = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/NICE_inverse_correct_boundary/rampup_5-15s","r")
# entry = imas.DBEntry("imas:hdf5?path=/work/imas/shared/imasdb/ITER/3/105084/1","r")



entry_machine_description = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/3","r") # TODO: take geometry information from a general source? Machine description database?

# Import IDSs that contain equilibrium information
eq = entry.get("equilibrium")                                          # get equilibirum data at requested time


# Import IDSs that contain geometry information
pf_active = entry.get("pf_active")                                     # get pf_active data (in principle also only geometry data of coils that is taken from this in this script)
wall = entry_machine_description.get("wall", lazy=True)                           # taken from other entry, since inverse solutions do not contain wall, iron_core, pf_passive


# Time data from magnetic equilibrium 
time = eq.time
number_of_time_slices = len(time)


# Coil parameters
number_of_coils = pf_active.coil.size
print("Number of coils = ", number_of_coils)


# The last two coils, Upper and Lower vertical Stabilization (VSL and VSU) coils, do not have rectangular shape, but the outline must be taken
VSU_outline_r = pf_active.coil[12].element[0].geometry.outline.r
VSU_outline_z = pf_active.coil[12].element[0].geometry.outline.z

VSL_outline_r = pf_active.coil[13].element[0].geometry.outline.r
VSL_outline_z = pf_active.coil[13].element[0].geometry.outline.z



# Compute global psi_min and psi_max across all time slices
psi_min = float("inf")
psi_max = float("-inf")



# When equilibrium data is under profiles_2d field: e.g. for DINA results
if equilibrium_data_stored == 'profiles_2d':
    for time_slice in eq.time_slice:
        psi_values = time_slice.profiles_2d[0].psi
        psi_min = min(psi_min, psi_values.min())
        psi_max = max(psi_max, psi_values.max())


# when equilibrium data is under ggd field: e.g. for NICE inverse results
if equilibrium_data_stored == 'ggd':
    for time_slice in eq.time_slice:
        psi_values = time_slice.ggd[0].psi[0].values
        psi_min = min(psi_min, psi_values.min())
        psi_max = max(psi_max, psi_values.max())

    


# The contour values of psi now remain constant throughout time for the slider plot
number_of_contours = 100                                    # Here specify the amount of psi iso contours that are visualised in the plot
contour_values_psi = np.linspace(psi_min, psi_max, number_of_contours) 




# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 10))
plt.subplots_adjust(bottom=0.3)  # Adjust for slider space


# Add the rectangular coil outlines to the plot
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
        coil_width,                          # Width
        coil_height,                         # Height
        edgecolor=(0, 0, 0, 0.8),               # Outline color
        facecolor="none",               # Transparent fill
        linewidth=2,                    # Outline thickness
    )
    ax.add_patch(rect)



# Plotting VS coils

ax.plot(VSU_outline_r, VSU_outline_z, 'o', markersize=0.6, color="black")
ax.plot(VSL_outline_r, VSL_outline_z, 'o', markersize=0.6, color="black")



# Plotting ITER inner and outer wall outline

r_wall_inner = wall.description_2d[0].vessel.unit[0].annular.centreline.r
z_wall_inner = wall.description_2d[0].vessel.unit[0].annular.centreline.z

r_wall_outer = wall.description_2d[0].vessel.unit[1].annular.centreline.r
z_wall_outer = wall.description_2d[0].vessel.unit[1].annular.centreline.z

ax.plot(r_wall_inner, z_wall_inner, color="black", label="Inner wall contour", alpha=0.8)
ax.plot(r_wall_outer, z_wall_outer, color="black", label="Outer wall contour", alpha=0.8)




# Plotting ITER divertor contour

r_limiter = wall.description_2d[0].limiter.unit[0].outline.r
z_limiter = wall.description_2d[0].limiter.unit[0].outline.z

r_divertor = wall.description_2d[0].limiter.unit[1].outline.r
z_divertor = wall.description_2d[0].limiter.unit[1].outline.z

ax.plot(r_limiter, z_limiter, color="black", label="Limiter contour", alpha=0.8)
ax.plot(r_divertor, z_divertor, color="black", label="Divertor contour", alpha=0.8)


# Initial plot data for initial time slice 
boundary_line, = ax.plot([], [], '.', color="red", markersize=2, label="Plasma boundary")
contour_plot = None
psi_plot = None

# Create text objects once
psi_axis_text = ax.text(
    -0.25, 0.8, "", ha="center", va="top",
    transform=ax.transAxes, fontsize=12, fontweight="bold"
)

psi_boundary_text = ax.text(
    -0.25, 0.7, "", ha="center", va="top",
    transform=ax.transAxes, fontsize=12, fontweight="bold"
)

# Formatting
ax.set_aspect('equal')
ax.set_xlabel("R (m)")
ax.set_ylabel("Z (m)")
ax.set_xlim(0, 13)
ax.set_ylim(-10, 10)




# Add slider
ax_slider = plt.axes([0.2, -0.008, 0.6, 0.03], facecolor="lightgrey")
time_slider = Slider(ax_slider, "Time", 0, number_of_time_slices - 1, valinit=0, valstep=1)


# Add a colorbar for psi scatter plot
norm = matplotlib.colors.Normalize(vmin=psi_min, vmax=psi_max)
colorbar = plt.colorbar(
    matplotlib.cm.ScalarMappable(norm=norm, cmap="viridis_r"),
    ax=ax,
    orientation="vertical",
    label=r"Poloidal Magnetic Flux ($\psi$)"
)


# Handle key press events for slider control
fig.canvas.mpl_connect('key_press_event', on_key)


def update_plot(i):
    """Update the plot for a given time index i."""
    global contour_plot, psi_plot
    time_idx = int(i)
    


    ## Get data for the selected time slice

    # Plasma boundary data
    plasma_boundary_r = eq.time_slice[time_idx].boundary.outline.r
    plasma_boundary_z = eq.time_slice[time_idx].boundary.outline.z
    



    # Extract poloidal magnetic flux data 

    # When equilibrium data is under profiles_2d field: e.g. for DINA results
    if equilibrium_data_stored == 'profiles_2d':
        r = eq.time_slice[time_idx].profiles_2d[0].r
        z = eq.time_slice[time_idx].profiles_2d[0].z
        psi = eq.time_slice[time_idx].profiles_2d[0].psi

    # when equilibrium data is under ggd field: e.g. for NICE inverse results
    if equilibrium_data_stored == 'ggd':
        r = eq.time_slice[time_idx].ggd[0].r[0].values
        z = eq.time_slice[time_idx].ggd[0].z[0].values
        psi = eq.time_slice[time_idx].ggd[0].psi[0].values


    # Interesting value of the poloidal magnetic flux (psi) are on the magnetic axis and plasma boundary - these will be displayed below the plots for the input and output data

    psi_axis = eq.time_slice[time_idx].global_quantities.psi_axis.value
    psi_boundary = eq.time_slice[time_idx].global_quantities.psi_boundary.value    
   

    # Update plasma boundary
    boundary_line.set_data(plasma_boundary_r, plasma_boundary_z)
    
    # Remove previous contour and scatter plot
    if contour_plot is not None:
        for coll in contour_plot.collections:
            coll.remove()

    if psi_plot is not None:
        psi_plot.remove()


    # Redraw contour and scatter plot
    contour_plot = ax.tricontour(r.flatten(), z.flatten(), psi.flatten(), levels=contour_values_psi, 
                                 linewidths=1, cmap="viridis_r")
    opaqueness_of_psi_scatter = 0.01
    psi_plot = ax.scatter(r, z, c=psi, cmap="viridis_r", alpha=opaqueness_of_psi_scatter)


    # Update title with time
    ax.set_title(f"Magnetic equilibrium poloidal flux contours at time = {time[time_idx]:.4f}s")


    # Formatting
    ax.set_aspect('equal')
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_xlim(0, 13)
    ax.set_ylim(-10, 10)
    
    psi_axis_text.set_text(rf"$ψ_{{axis}}$={psi_axis:.3e} Wb")
    psi_boundary_text.set_text(rf"$ψ_{{boundary}}$={psi_boundary:.3e} Wb")
    


    # Refresh plot
    fig.canvas.draw_idle()

# Connect slider to update function
time_slider.on_changed(update_plot)

# Initial plot update
update_plot(0)

plt.tight_layout()
plt.show()