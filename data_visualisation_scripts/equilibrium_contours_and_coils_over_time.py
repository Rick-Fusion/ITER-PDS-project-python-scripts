import imas
import imaspy
import imaspy.util
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from matplotlib.widgets import Slider



# function for controlling slider with left and right keys

def on_key(event):
    """Handle key presses to adjust the slider value."""
    if event.key == 'right':
        # Increase slider value by 1 (move right)
        time_slider.set_val(min(time_slider.val + 1, time_slider.valmax))
    elif event.key == 'left':
        # Decrease slider value by 1 (move left)
        time_slider.set_val(max(time_slider.val - 1, time_slider.valmin))



# This python script is dedicated to plotting contours of the poloidal magnetic flux (psi) in a 2D plane from the equilibrium IDS over the duration of equilibrium data 
# plus PF and CS coil cross-sections



matplotlib.use("TkAgg")

entry = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_19/iteration_1/post_inverse","r")

# Import IDSs
eq = entry.get("equilibrium")         # get equilibirum data at requested time
pf_active = entry.get("pf_active")                                                          # get pf_active data


# Time data
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

'''
# only for DINA data (run 1 & 3 & 4)
for time_slice in eq.time_slice:
    psi_values = time_slice.profiles_2d[0].psi
    psi_min = min(psi_min, psi_values.min())
    psi_max = max(psi_max, psi_values.max())
'''
test=0
# For Inverse results
for time_slice in eq.time_slice:
    test+=1
    print(test)
    psi_values = time_slice.ggd[0].psi[0].values
    psi_min = min(psi_min, psi_values.min())
    psi_max = max(psi_max, psi_values.max())

    
# The contour values of psi now remain constant throughout time for the slider plot
contour_values_psi = np.arange(psi_min, psi_max, 10) 



# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 10))
plt.subplots_adjust(bottom=0.3)  # Adjust for slider space


# Add coils as rectangles
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
        edgecolor="black",               # Outline color
        facecolor="none",               # Transparent fill
        linewidth=2                     # Outline thickness
    )
    ax.add_patch(rect)


# Plotting VS coils

ax.plot(VSU_outline_r, VSU_outline_z, 'o', markersize=0.6, color="black")
ax.plot(VSL_outline_r, VSL_outline_z, 'o', markersize=0.6, color="black")



# Initial plot data for time slice 0
boundary_line, = ax.plot([], [], '-', color="red", label="Plasma boundary")
contour_plot = None
psi_plot = None

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
    
    # Get data for the selected time slice

    # Plasma boundary data
    plasma_boundary_r = eq.time_slice[time_idx].boundary.outline.r
    plasma_boundary_z = eq.time_slice[time_idx].boundary.outline.z
    
    # Poloidal magnetic flux data

    '''
    # only for DINA data (run 1 & 3 & 4)
    r = eq.time_slice[time_idx].profiles_2d[0].r
    z = eq.time_slice[time_idx].profiles_2d[0].z
    psi = eq.time_slice[time_idx].profiles_2d[0].psi
    '''
    
    # For Inverse results
    r = eq.time_slice[time_idx].ggd[0].r[0].values
    z = eq.time_slice[time_idx].ggd[0].z[0].values
    psi = eq.time_slice[time_idx].ggd[0].psi[0].values
        

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
    psi_plot = ax.scatter(r, z, c=psi, cmap="viridis_r", alpha=0.01)


    # Update title with time
    ax.set_title(f"Equilibrium poloidal flux contours at time = {time[time_idx]:.2f}s + PF Coil Cross-Sections")


    # Formatting
    ax.set_aspect('equal')
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_xlim(0, 13)
    ax.set_ylim(-10, 10)
    
    # Refresh plot
    fig.canvas.draw_idle()

# Connect slider to update function
time_slider.on_changed(update_plot)

# Initial plot update
update_plot(0)

plt.tight_layout()
plt.show()