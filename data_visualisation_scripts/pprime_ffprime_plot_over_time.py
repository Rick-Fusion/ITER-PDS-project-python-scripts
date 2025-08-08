import imas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
matplotlib.use("TkAgg")

"""The aim of this plotting script is to visualise the pprime and ffprime profiles over time for multiple iterations resulting from the torax ouput or from any point in the NICE inverse + TORAX coupling scheme"""


### NOTE: MAKE SURE THAT ALL EQUILIBRIUM IDS's ARE CONSISTENT IN THE TIME FOR EACH TIME SLICE, SINCE THEY ARE COMPARED IN THE SAME FIGURE FOR EACH TIME SLICE


def on_key(event):
    """Handle key presses to adjust the slider value."""
    if event.key == 'right':
        # Increase slider value by 1 (move right)
        time_slider.set_val(min(time_slider.val + 1, time_slider.valmax))
    elif event.key == 'left':
        # Decrease slider value by 1 (move left)
        time_slider.set_val(max(time_slider.val - 1, time_slider.valmin))


#open the Data Entry
#entry_0 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/5","r")
entry_0 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_2/iteration_1/post_inverse", "r")
entry_1 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_2/iteration_1/post_torax", "r")
#entry_2 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_2/iteration_1/post_torax", "r")
#entry_3 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_2/iteration_1/post_torax", "r")
#entry_4 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_2/iteration_1/post_torax", "r")
#entry_5 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_2/iteration_1/post_torax", "r")

# get the full equilibrium IDS for all the iteration for 1 and higher
eq_1 = entry_1.get("equilibrium", lazy=True)
#eq_2 = entry_2.get("equilibrium", lazy=True)
#eq_3 = entry_3.get("equilibrium", lazy=True)
#eq_4 = entry_4.get("equilibrium", lazy=True)
#eq_5 = entry_5.get("equilibrium", lazy=True)

# Time data
time_1 = eq_1.time
#time_2 = eq_2.time
#time_3 = eq_3.time
#time_4 = eq_4.time
#time_5 = eq_5.time


number_of_time_slices = len(time_1)


# For entry_0 (the initial data), the times do not match, since only a certain interval is taken
# Therefore create the eq_0 for which time times do match with the other eq's
entry_0_mem = imas.DBEntry("imas:memory?path=/", "w")

for time in time_1:
    eq_0_slice = entry_0.get_slice("equilibrium", time, imas.ids_defs.CLOSEST_INTERP)
    entry_0_mem.put_slice(eq_0_slice)

eq_0 = entry_0_mem.get("equilibrium", lazy=True)

time_0 = eq_0.time

print(f"## Time slices in eq_0 = {len(time_0)}##")
print(f"## Time slices in eq_1 = {len(time_1)}##")
#print(f"## Time slices in eq_2 = {len(time_2)}##")
#print(f"## Time slices in eq_3 = {len(time_3)}##")
#print(f"## Time slices in eq_4 = {len(time_4)}##")
#print(f"## Time slices in eq_5 = {len(time_5)}##")



# Create a figure and axis
fig, (ax1 , ax2) = plt.subplots(1,2, figsize=(12, 10))
plt.subplots_adjust(bottom=0.25)  # Adjust for slider space

# Initial plot data for time slice 0
pprime_plot_0, = ax1.plot([], [], '-', color="green", label="pprime profile it=0")
ffprime_plot_0, = ax2.plot([], [], '-', color="green", label="ffrime profile it=0")

pprime_plot_1, = ax1.plot([], [], '-', color="red", label="pprime profile it=1")
ffprime_plot_1, = ax2.plot([], [], '-', color="red", label="ffrime profile it=1")

pprime_plot_2, = ax1.plot([], [], '-', color="blue", label="pprime profile it=2")
ffprime_plot_2, = ax2.plot([], [], '-', color="blue", label="ffrime profile it=2")

pprime_plot_3, = ax1.plot([], [], '-', color="magenta", label="pprime profile it=3")
ffprime_plot_3, = ax2.plot([], [], '-', color="magenta", label="ffrime profile it=3")

pprime_plot_4, = ax1.plot([], [], '-', color='yellow', label="pprime profile it=4")
ffprime_plot_4, = ax2.plot([], [], '-', color='yellow', label="ffrime profile it=4")

pprime_plot_5, = ax1.plot([], [], '-', color='black', label="pprime profile it=5")
ffprime_plot_5, = ax2.plot([], [], '-', color='black', label="ffrime profile it=5")

# Add slider
ax_slider = plt.axes([0.2, -0.005, 0.6, 0.03], facecolor="lightgrey")
time_slider = Slider(ax_slider, "Time", 0, number_of_time_slices - 1, valinit=0, valstep=1)

# Handle key press events for slider control
fig.canvas.mpl_connect('key_press_event', on_key)

def update_plot(i):
    """Update the plot for a given time index i."""
    time_idx = int(i)
    
    # Get data for the selected time slice

    # Profiles data

    pprime_0 = eq_0.time_slice[time_idx].profiles_1d.dpressure_dpsi
    ffprime_0 = eq_0.time_slice[time_idx].profiles_1d.f_df_dpsi
    psi_0 = eq_0.time_slice[time_idx].profiles_1d.psi
    psi_norm_0 = ( psi_0 - psi_0[0] ) / ( psi_0[-1] - psi_0[0] )
    rho_tor_norm_0 = eq_0.time_slice[time_idx].profiles_1d.rho_tor_norm
    pressure_0 = eq_0.time_slice[time_idx].profiles_1d.pressure


    pprime_1 = eq_1.time_slice[time_idx].profiles_1d.dpressure_dpsi
    ffprime_1 = eq_1.time_slice[time_idx].profiles_1d.f_df_dpsi
    psi_1 = eq_1.time_slice[time_idx].profiles_1d.psi
    psi_norm_1 = ( psi_1 - psi_1[0] ) / ( psi_1[-1] - psi_1[0] )
    rho_tor_norm_1 = eq_1.time_slice[time_idx].profiles_1d.rho_tor_norm
    pressure_1 = eq_1.time_slice[time_idx].profiles_1d.pressure
    
    """
    pprime_2 = eq_2.time_slice[time_idx].profiles_1d.dpressure_dpsi
    ffprime_2 = eq_2.time_slice[time_idx].profiles_1d.f_df_dpsi
    psi_2 = eq_2.time_slice[time_idx].profiles_1d.psi
    psi_norm_2 = ( psi_2 - psi_2[0] ) / ( psi_2[-1] - psi_2[0] )
    rho_tor_norm_2 = eq_2.time_slice[time_idx].profiles_1d.rho_tor_norm
    pressure_2 = eq_2.time_slice[time_idx].profiles_1d.pressure
   

    pprime_3 = eq_3.time_slice[time_idx].profiles_1d.dpressure_dpsi
    ffprime_3 = eq_3.time_slice[time_idx].profiles_1d.f_df_dpsi
    psi_3 = eq_3.time_slice[time_idx].profiles_1d.psi
    psi_norm_3 = ( psi_3 - psi_3[0] ) / ( psi_3[-1] - psi_3[0] )
    rho_tor_norm_3 = eq_3.time_slice[time_idx].profiles_1d.rho_tor_norm
    pressure_3 = eq_3.time_slice[time_idx].profiles_1d.pressure

    pprime_4 = eq_4.time_slice[time_idx].profiles_1d.dpressure_dpsi
    ffprime_4 = eq_4.time_slice[time_idx].profiles_1d.f_df_dpsi
    psi_4 = eq_4.time_slice[time_idx].profiles_1d.psi
    psi_norm_4 = ( psi_4 - psi_4[0] ) / ( psi_4[-1] - psi_4[0] )
    rho_tor_norm_4 = eq_4.time_slice[time_idx].profiles_1d.rho_tor_norm
    pressure_4 = eq_4.time_slice[time_idx].profiles_1d.pressure

    pprime_5 = eq_5.time_slice[time_idx].profiles_1d.dpressure_dpsi
    ffprime_5 = eq_5.time_slice[time_idx].profiles_1d.f_df_dpsi
    psi_5 = eq_5.time_slice[time_idx].profiles_1d.psi
    psi_norm_5 = ( psi_5 - psi_5[0] ) / ( psi_5[-1] - psi_5[0] )
    rho_tor_norm_5 = eq_5.time_slice[time_idx].profiles_1d.rho_tor_norm
    pressure_5 = eq_5.time_slice[time_idx].profiles_1d.pressure
    """

    # Update profiles
    pprime_plot_0.set_data(psi_norm_0, pprime_0)
    ffprime_plot_0.set_data(psi_norm_0, ffprime_0)
    
    pprime_plot_1.set_data(psi_norm_1, pprime_1)
    ffprime_plot_1.set_data(psi_norm_1, ffprime_1)
    """
    pprime_plot_2.set_data(psi_norm_2, pprime_2)
    ffprime_plot_2.set_data(psi_norm_2, ffprime_2)

    pprime_plot_3.set_data(psi_norm_3, pprime_3)
    ffprime_plot_3.set_data(psi_norm_3, ffprime_3)

    pprime_plot_4.set_data(psi_norm_4, pprime_4)
    ffprime_plot_4.set_data(psi_norm_4, ffprime_4)

    pprime_plot_5.set_data(psi_norm_5, pprime_5)
    ffprime_plot_5.set_data(psi_norm_5, ffprime_5)
    """
    # Update title with time
    ax1.set_title(f"p' profile at time = {time_1[time_idx]:.2f}s")     # optional: from torax output
    ax2.set_title(f"ff' profile at time = {time_1[time_idx]:.2f}s")    # optional: from torax output


    # Formatting
    #ax1.set_aspect('equal')
    ax1.set_xlabel(r"$\psi_{norm}$")
    ax1.set_ylabel(r"p' (Pa $Wb^{-1}$)")
    #ax1.set_xlim(0, 1)
    #ax1.set_ylim(-1, 1)
    ax1.relim()
    ax1.autoscale_view(True, True, True)
    ax1.legend(fontsize=7, loc="lower right", frameon=True)
    
    # Formatting
    #ax2.set_aspect('equal')
    ax2.set_xlabel(r"$\psi_{norm}$")
    ax2.set_ylabel(r"ff' ")
    #ax2.set_xlim(0, 1)
    #ax2.set_ylim(-1, 1)
    ax2.relim()
    ax2.autoscale_view(True, True, True)
    ax2.legend(fontsize=7, loc="upper right", frameon=True)

    # Refresh plot
    fig.canvas.draw_idle()

# Connect slider to update function
time_slider.on_changed(update_plot)

# Initial plot update
update_plot(0)

plt.tight_layout()
plt.show()

entry_0_mem.close()
entry_0.close()
entry_1.close()
#entry_2.close()
#entry_3.close()
#entry_4.close()
#entry_5.close()