import imas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
matplotlib.use("TkAgg")

"""The aim of this plotting script is to visualise the equilibrium plasma boundary outline over time resulting from the torax ouput in the NICE inverse + TORAX coupling scheme"""

def on_key(event):
    """Handle key presses to adjust the slider value."""
    if event.key == 'right':
        # Increase slider value by 1 (move right)
        time_slider.set_val(min(time_slider.val + 1, time_slider.valmax))
    elif event.key == 'left':
        # Decrease slider value by 1 (move left)
        time_slider.set_val(max(time_slider.val - 1, time_slider.valmin))


# open the Data Entry
entry = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_8/post_torax","r")

# get the full equilibrium IDS
eq = entry.get("equilibrium")

# Time data
time = eq.time
number_of_time_slices = len(time)


# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 10))
plt.subplots_adjust(bottom=0.3)  # Adjust for slider space

# Initial plot data for time slice 0
boundary_line, = ax.plot([], [], '-', color="red", label="Plasma boundary")


# Add slider
ax_slider = plt.axes([0.2, -0.008, 0.6, 0.03], facecolor="lightgrey")
time_slider = Slider(ax_slider, "Time", 0, number_of_time_slices - 1, valinit=0, valstep=1)

# Handle key press events for slider control
fig.canvas.mpl_connect('key_press_event', on_key)

def update_plot(i):
    """Update the plot for a given time index i."""
    time_idx = int(i)
    
    # Get data for the selected time slice

    # Plasma boundary data
    plasma_boundary_r = eq.time_slice[time_idx].boundary.outline.r
    plasma_boundary_z = eq.time_slice[time_idx].boundary.outline.z
    
   
    # Update plasma boundary
    boundary_line.set_data(plasma_boundary_r, plasma_boundary_z)

    # Update title with time
    ax.set_title(f"Equilibrium plasma boudnary from torax output at time = {time[time_idx]:.2f}s")


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