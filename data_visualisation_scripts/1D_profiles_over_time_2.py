import imas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
matplotlib.use("TkAgg")

""" 
The aim of this plotting script is to visualise 1D profiles in general over time for multiple iterations resulting from the torax ouput or from any point in the NICE inverse + TORAX coupling scheme 

However, this plotting script is mainly used to 1D profiles between input and output of TORAX to check consistency. You can plot these quantities as well over the whole workflow as is desired
If you would like to see p' and ff' convergence between iterations, you should look at 1D_profiles_over_time_1


The 1D profiles that are currently plotted in this script include: 
- psi_norm vs rho_tor_norm    -   poloidal magnetic flux normalised
- j_phi vs rho_tor_norm    -   current density profile
- phi vs rho_tor_norm   -    toroidal magnetic flux
- gm1, gm2, gm3, gm7 vs rho_tor_norm    -    geometric parameters
- dpsi_drho_tor vs rho_tor_norm
- dvolume_dpsi vs rho_tor_norm
- q vs rho_tor_norm  -   safety factor


"""

### NOTE: MAKE SURE THAT ALL EQUILIBRIUM IDS's ARE CONSISTENT IN THE TIME FOR EACH TIME SLICE, SINCE THEY ARE COMPARED IN THE SAME FIGURE FOR EACH TIME SLICE


def on_key(event):
    """Handle key presses to adjust the slider value."""
    if event.key == 'right':
        # Increase slider value by 1 (move right)
        time_slider.set_val(min(time_slider.val + 1, time_slider.valmax))
    elif event.key == 'left':
        # Decrease slider value by 1 (move left)
        time_slider.set_val(max(time_slider.val - 1, time_slider.valmin))



###MAYBE TO-DO Potential easier implementation for selecting data entries below
# for path in [
#     'path1',
#     'path1',
#     # 'path1',
#     # 'path1',
#     # 'path1',
#     'path1',
# ]:
# if use_entry_0:
#     entry_0 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/initial_data","r")
#     eq_0 = entry_0.get("equilibrium", lazy=True)
#     time_0 = eq_0.time
#     print(f"## Time slices in eq_0 = {len(time_0)}##")
#     number_of_time_slices = len(time_0)
#     time = time_0



# Here you can specify which entry to use and then then all the lines below will be selected accordingly, so you don't have to manually comment lines out manually each time :)
use_entry_0 = False
use_entry_1 = True      
use_entry_2 = True
use_entry_3 = False
use_entry_4 = False
use_entry_5 = False
use_entry_6 = False
use_entry_7 = False
use_entry_8 = False
use_entry_9 = False
use_entry_10 = False


#open the Data Entry and get equilibrium
if use_entry_0:
    entry_0 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7.2/initial_data","r")
    eq_0 = entry_0.get("equilibrium", lazy=True)
    time_0 = eq_0.time
    print(f"## Time slices in eq_0 = {len(time_0)}##")

if use_entry_1:
    entry_1 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_10/iteration_1/post_inverse","r")
    eq_1 = entry_1.get("equilibrium", lazy=True)
    time_1 = eq_1.time
    print(f"## Time slices in eq_1 = {len(time_1)}##")
if use_entry_2:
    entry_2 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_10/iteration_1/post_torax","r")
    eq_2 = entry_2.get("equilibrium", lazy=True)
    time_2 = eq_2.time
    print(f"## Time slices in eq_2 = {len(time_2)}##")

if use_entry_3:
    entry_3 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_3/post_torax","r")
    eq_3 = entry_3.get("equilibrium", lazy=True)
    time_3 = eq_3.time
    print(f"## Time slices in eq_3 = {len(time_3)}##")
if use_entry_4:
    entry_4 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_4/post_torax","r")
    eq_4 = entry_4.get("equilibrium", lazy=True)
    time_4 = eq_4.time
    print(f"## Time slices in eq_4 = {len(time_4)}##")
if use_entry_5:
    entry_5 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_5/post_torax","r")
    eq_5 = entry_5.get("equilibrium", lazy=True)
    time_5 = eq_5.time
    print(f"## Time slices in eq_5 = {len(time_5)}##")
if use_entry_6:
    entry_6 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_6/post_torax","r")
    eq_6 = entry_6.get("equilibrium", lazy=True)
    time_6 = eq_6.time
    print(f"## Time slices in eq_6 = {len(time_6)}##")
if use_entry_7:
    entry_7 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_7/post_torax","r")
    eq_7 = entry_7.get("equilibrium", lazy=True)
    time_7 = eq_7.time
    print(f"## Time slices in eq_7 = {len(time_7)}##")
if use_entry_8:
    entry_8 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_8/post_torax","r")
    eq_8 = entry_8.get("equilibrium", lazy=True)
    time_8 = eq_8.time
    print(f"## Time slices in eq_8 = {len(time_8)}##")
if use_entry_9:
    entry_9 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_9/post_torax","r")
    eq_9 = entry_9.get("equilibrium", lazy=True)
    time_9 = eq_9.time
    print(f"## Time slices in eq_9 = {len(time_9)}##")
if use_entry_10:
    entry_10 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_7/iteration_10/post_torax","r")
    eq_10 = entry_10.get("equilibrium", lazy=True)
    time_10 = eq_10.time
    print(f"## Time slices in eq_10 = {len(time_10)}##")


if use_entry_0:
    number_of_time_slices = len(time_0)
    time = time_0
elif use_entry_1:
    number_of_time_slices = len(time_1)
    time = time_1
elif use_entry_2:
    number_of_time_slices = len(time_2)
    time = time_2
elif use_entry_3:
    number_of_time_slices = len(time_3)
    time = time_3
elif use_entry_4:
    number_of_time_slices = len(time_4)
    time = time_4
elif use_entry_5:
    number_of_time_slices = len(time_5)
    time = time_5
elif use_entry_6:
    number_of_time_slices = len(time_6)
    time = time_6
elif use_entry_7:
    number_of_time_slices = len(time_7)
    time = time_7
elif use_entry_8:
    number_of_time_slices = len(time_8)
    time = time_8
elif use_entry_9:
    number_of_time_slices = len(time_9)
    time = time_9
elif use_entry_10:
    number_of_time_slices = len(time_10)
    time = time_10

"""
# For entry_0 (the initial data), the times do not match, since only a certain interval is taken
# Therefore create the eq_0 for which time times do match with the other eq's
entry_0_mem = imas.DBEntry("imas:memory?path=/", "w")

for time in time_1:
    eq_0_slice = entry_0.get_slice("equilibrium", time, imas.ids_defs.CLOSEST_INTERP)
    entry_0_mem.put_slice(eq_0_slice)

eq_0 = entry_0_mem.get("equilibrium", lazy=True)

time_0 = eq_0.time
"""




# Create a figure and axis
fig, ((ax1, ax2, ax3, ax4, ax5),(ax6, ax7, ax8, ax9, ax10)) = plt.subplots(2,5, figsize=(20, 20), gridspec_kw={'wspace': 0.5, 'hspace': 0.35})
#plt.subplots_adjust(bottom=0.25)  # Adjust for slider space


# Initial plot data for time slice 0

if use_entry_0:
    psi_norm_plot_0, = ax1.plot([], [], '-', color="green", label="psi_norm profile it=0")

if use_entry_1:
    psi_norm_plot_1, = ax1.plot([], [], '-', color="red", label="post_inverse")
    j_phi_plot_1, = ax2.plot([], [], '-', color="red", label="post_inverse")
    phi_plot_1, = ax3.plot([], [], '-', color="red", label="post_inverse")
    q_plot_1, = ax4.plot([], [], '-', color="red", label="post_inverse")
    dpsi_drho_tor_plot_1, = ax5.plot([], [], '-', color="red", label="post_inverse")
    gm1_plot_1, = ax6.plot([], [], '-', color="red", label="post_inverse")
    gm2_plot_1, = ax7.plot([], [], '-', color="red", label="post_inverse")
    gm3_plot_1, = ax8.plot([], [], '-', color="red", label="post_inverse")
    gm7_plot_1, = ax9.plot([], [], '-', color="red", label="post_inverse")
    dvolume_dpsi_plot_1, = ax10.plot([], [], '-', color="red", label="post_inverse")
    

if use_entry_2:
    psi_norm_plot_2, = ax1.plot([], [], '-', color="blue", label="post_torax")
    j_phi_plot_2, = ax2.plot([], [], '-', color="blue", label="post_torax")
    phi_plot_2, = ax3.plot([], [], '-', color="blue", label="post_torax")
    q_plot_2, = ax4.plot([], [], '-', color="blue", label="post_torax")
    dpsi_drho_tor_plot_2, = ax5.plot([], [], '-', color="blue", label="post_torax")
    gm1_plot_2, = ax6.plot([], [], '-', color="blue", label="post_torax")
    gm2_plot_2, = ax7.plot([], [], '-', color="blue", label="post_torax")
    gm3_plot_2, = ax8.plot([], [], '-', color="blue", label="post_torax")
    gm7_plot_2, = ax9.plot([], [], '-', color="blue", label="post_torax")
    dvolume_dpsi_plot_2, = ax10.plot([], [], '-', color="blue", label="post_torax")



if use_entry_3:
    pprime_plot_3, = ax1.plot([], [], '-', color="magenta", label="pprime profile it=3")
    ffprime_plot_3, = ax2.plot([], [], '-', color="magenta", label="ffrime profile it=3")
    pressure_plot_3, = ax3.plot([], [], '-', color="magenta", label="pressure profile it=3")
    psi_norm_plot_3, = ax4.plot([], [], '-', color="magenta", label="psi_norm profile it=3")
 
if use_entry_4:
    pprime_plot_4, = ax1.plot([], [], '-', color='yellow', label="pprime profile it=4")
    ffprime_plot_4, = ax2.plot([], [], '-', color='yellow', label="ffrime profile it=4")
    pressure_plot_4, = ax3.plot([], [], '-', color="yellow", label="pressure profile it=4")
    psi_norm_plot_4, = ax4.plot([], [], '-', color="yellow", label="psi_norm profile it=4")

if use_entry_5:
    pprime_plot_5, = ax1.plot([], [], '-', color='black', label="pprime profile it=5")
    ffprime_plot_5, = ax2.plot([], [], '-', color='black', label="ffrime profile it=5")
    pressure_plot_5, = ax3.plot([], [], '-', color="black", label="pressure profile it=5")
    psi_norm_plot_5, = ax4.plot([], [], '-', color="black", label="psi_norm profile it=5")

if use_entry_6:
    pprime_plot_6, = ax1.plot([], [], '-',  label="pprime profile it=6")
    ffprime_plot_6, = ax2.plot([], [], '-',  label="ffrime profile it=6")
    pressure_plot_6, = ax3.plot([], [], '-',  label="pressure profile it=6")
    psi_norm_plot_6, = ax4.plot([], [], '-', label="psi_norm profile it=6")

if use_entry_7:
    pprime_plot_7, = ax1.plot([], [], '-', label="pprime profile it=7")
    ffprime_plot_7, = ax2.plot([], [], '-',  label="ffrime profile it=7")
    pressure_plot_7, = ax3.plot([], [], '-',  label="pressure profile it=7")
    psi_norm_plot_7, = ax4.plot([], [], '-',  label="psi_norm profile it=7")
 
if use_entry_8:
    pprime_plot_8, = ax1.plot([], [], '-',  label="pprime profile it=8")
    ffprime_plot_8, = ax2.plot([], [], '-',  label="ffrime profile it=8")
    pressure_plot_8, = ax3.plot([], [], '-',  label="pressure profile it=8")
    psi_norm_plot_8, = ax4.plot([], [], '-', label="psi_norm profile it=8")

if use_entry_9:
    pprime_plot_9, = ax1.plot([], [], '-',  label="pprime profile it=9")
    ffprime_plot_9, = ax2.plot([], [], '-',  label="ffrime profile it=9")
    pressure_plot_9, = ax3.plot([], [], '-',  label="pressure profile it=9")
    psi_norm_plot_9, = ax4.plot([], [], '-',  label="psi_norm profile it=9")

if use_entry_10:
    pprime_plot_10, = ax1.plot([], [], '-', label="pprime profile it=10")
    ffprime_plot_10, = ax2.plot([], [], '-', label="ffrime profile it=10")
    pressure_plot_10, = ax3.plot([], [], '-', label="pressure profile it=10")
    psi_norm_plot_10, = ax4.plot([], [], '-', label="psi_norm profile it=10")



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
    if use_entry_0:
        pprime_0 = eq_0.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_0 = eq_0.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_0 = eq_0.time_slice[time_idx].profiles_1d.psi
        psi_norm_0 = ( psi_0 - psi_0[0] ) / ( psi_0[-1] - psi_0[0] )
        rho_tor_norm_0 = eq_0.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_0 = eq_0.time_slice[time_idx].profiles_1d.pressure

    if use_entry_1:
        pprime_1 = eq_1.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_1 = eq_1.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_1 = eq_1.time_slice[time_idx].profiles_1d.psi
        psi_norm_1 = ( psi_1 - psi_1[0] ) / ( psi_1[-1] - psi_1[0] )
        rho_tor_norm_1 = eq_1.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_1 = eq_1.time_slice[time_idx].profiles_1d.pressure
        j_phi_1 = eq_1.time_slice[time_idx].profiles_1d.j_phi
        phi_1 = eq_1.time_slice[time_idx].profiles_1d.phi
        q_1 = eq_1.time_slice[time_idx].profiles_1d.q
        dpsi_drho_tor_1 = eq_1.time_slice[time_idx].profiles_1d.dpsi_drho_tor
        gm1_1 = eq_1.time_slice[time_idx].profiles_1d.gm1
        gm2_1 = eq_1.time_slice[time_idx].profiles_1d.gm2
        gm3_1 = eq_1.time_slice[time_idx].profiles_1d.gm3
        gm7_1 = eq_1.time_slice[time_idx].profiles_1d.gm7
        dvolume_dpsi_1 = eq_1.time_slice[time_idx].profiles_1d.dvolume_dpsi

    if use_entry_2:
        pprime_2 = eq_2.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_2 = eq_2.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_2 = eq_2.time_slice[time_idx].profiles_1d.psi
        psi_norm_2 = ( psi_2 - psi_2[0] ) / ( psi_2[-1] - psi_2[0] )
        rho_tor_norm_2 = eq_2.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_2 = eq_2.time_slice[time_idx].profiles_1d.pressure
        j_phi_2 = eq_2.time_slice[time_idx].profiles_1d.j_phi
        phi_2 = eq_2.time_slice[time_idx].profiles_1d.phi
        q_2 = eq_2.time_slice[time_idx].profiles_1d.q
        dpsi_drho_tor_2 = eq_2.time_slice[time_idx].profiles_1d.dpsi_drho_tor
        gm1_2 = eq_2.time_slice[time_idx].profiles_1d.gm1
        gm2_2 = eq_2.time_slice[time_idx].profiles_1d.gm2
        gm3_2 = eq_2.time_slice[time_idx].profiles_1d.gm3
        gm7_2 = eq_2.time_slice[time_idx].profiles_1d.gm7
        dvolume_dpsi_2 = eq_2.time_slice[time_idx].profiles_1d.dvolume_dpsi

    if use_entry_3:
        pprime_3 = eq_3.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_3 = eq_3.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_3 = eq_3.time_slice[time_idx].profiles_1d.psi
        psi_norm_3 = ( psi_3 - psi_3[0] ) / ( psi_3[-1] - psi_3[0] )
        rho_tor_norm_3 = eq_3.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_3 = eq_3.time_slice[time_idx].profiles_1d.pressure

    if use_entry_4:
        pprime_4 = eq_4.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_4 = eq_4.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_4 = eq_4.time_slice[time_idx].profiles_1d.psi
        psi_norm_4 = ( psi_4 - psi_4[0] ) / ( psi_4[-1] - psi_4[0] )
        rho_tor_norm_4 = eq_4.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_4 = eq_4.time_slice[time_idx].profiles_1d.pressure

    if use_entry_5:
        pprime_5 = eq_5.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_5 = eq_5.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_5 = eq_5.time_slice[time_idx].profiles_1d.psi
        psi_norm_5 = ( psi_5 - psi_5[0] ) / ( psi_5[-1] - psi_5[0] )
        rho_tor_norm_5 = eq_5.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_5 = eq_5.time_slice[time_idx].profiles_1d.pressure

    if use_entry_6:
        pprime_6 = eq_6.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_6 = eq_6.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_6 = eq_6.time_slice[time_idx].profiles_1d.psi
        psi_norm_6 = ( psi_6 - psi_6[0] ) / ( psi_6[-1] - psi_6[0] )
        rho_tor_norm_6 = eq_6.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_6 = eq_6.time_slice[time_idx].profiles_1d.pressure

    if use_entry_7:
        pprime_7 = eq_7.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_7 = eq_7.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_7 = eq_7.time_slice[time_idx].profiles_1d.psi
        psi_norm_7 = ( psi_7 - psi_7[0] ) / ( psi_7[-1] - psi_7[0] )
        rho_tor_norm_7 = eq_7.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_7 = eq_7.time_slice[time_idx].profiles_1d.pressure

    if use_entry_8:
        pprime_8 = eq_8.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_8 = eq_8.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_8 = eq_8.time_slice[time_idx].profiles_1d.psi
        psi_norm_8 = ( psi_8 - psi_8[0] ) / ( psi_8[-1] - psi_8[0] )
        rho_tor_norm_8 = eq_8.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_8 = eq_8.time_slice[time_idx].profiles_1d.pressure

    if use_entry_9:
        pprime_9 = eq_9.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_9 = eq_9.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_9 = eq_9.time_slice[time_idx].profiles_1d.psi
        psi_norm_9 = ( psi_9 - psi_9[0] ) / ( psi_9[-1] - psi_9[0] )
        rho_tor_norm_9 = eq_9.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_9 = eq_9.time_slice[time_idx].profiles_1d.pressure

    if use_entry_10:
        pprime_10 = eq_10.time_slice[time_idx].profiles_1d.dpressure_dpsi
        ffprime_10 = eq_10.time_slice[time_idx].profiles_1d.f_df_dpsi
        psi_10 = eq_10.time_slice[time_idx].profiles_1d.psi
        psi_norm_10 = ( psi_10 - psi_10[0] ) / ( psi_10[-1] - psi_10[0] )
        rho_tor_norm_10 = eq_10.time_slice[time_idx].profiles_1d.rho_tor_norm
        pressure_10 = eq_10.time_slice[time_idx].profiles_1d.pressure
    
    

    # Update profiles
    if use_entry_0:
        psi_norm_plot_0.set_data(rho_tor_norm_0, psi_norm_0)
    

    if use_entry_1:
        psi_norm_plot_1.set_data(rho_tor_norm_1, psi_norm_1)
        j_phi_plot_1.set_data(rho_tor_norm_1, j_phi_1/1e6)              # in MA
        phi_plot_1.set_data(rho_tor_norm_1, phi_1)
        q_plot_1.set_data(rho_tor_norm_1, q_1)
        dpsi_drho_tor_plot_1.set_data(rho_tor_norm_1, dpsi_drho_tor_1)
        gm1_plot_1.set_data(rho_tor_norm_1, gm1_1)
        gm2_plot_1.set_data(rho_tor_norm_1, gm2_1)
        gm3_plot_1.set_data(rho_tor_norm_1, gm3_1)
        gm7_plot_1.set_data(rho_tor_norm_1, gm7_1)
        dvolume_dpsi_plot_1.set_data(rho_tor_norm_1, dvolume_dpsi_1)

    if use_entry_2:
        psi_norm_plot_2.set_data(rho_tor_norm_2, psi_norm_2)
        j_phi_plot_2.set_data(rho_tor_norm_2, j_phi_2/1e6)              # in MA
        phi_plot_2.set_data(rho_tor_norm_2, phi_2)
        q_plot_2.set_data(rho_tor_norm_2, q_2)
        dpsi_drho_tor_plot_2.set_data(rho_tor_norm_2, dpsi_drho_tor_2)
        gm1_plot_2.set_data(rho_tor_norm_2, gm1_2)
        gm2_plot_2.set_data(rho_tor_norm_2, gm2_2)
        gm3_plot_2.set_data(rho_tor_norm_2, gm3_2)
        gm7_plot_2.set_data(rho_tor_norm_2, gm7_2)
        dvolume_dpsi_plot_2.set_data(rho_tor_norm_2, dvolume_dpsi_2)


    if use_entry_3:
        pprime_plot_3.set_data(psi_norm_3, pprime_3)
        ffprime_plot_3.set_data(psi_norm_3, ffprime_3)
        pressure_plot_3.set_data(psi_norm_3, pressure_3)
        psi_norm_plot_3.set_data(rho_tor_norm_3, psi_norm_3)

    if use_entry_4:
        pprime_plot_4.set_data(psi_norm_4, pprime_4)
        ffprime_plot_4.set_data(psi_norm_4, ffprime_4)
        pressure_plot_4.set_data(psi_norm_4, pressure_4)
        psi_norm_plot_4.set_data(rho_tor_norm_4, psi_norm_4)

    if use_entry_5:
        pprime_plot_5.set_data(psi_norm_5, pprime_5)
        ffprime_plot_5.set_data(psi_norm_5, ffprime_5)
        pressure_plot_5.set_data(psi_norm_5, pressure_5)
        psi_norm_plot_5.set_data(rho_tor_norm_5, psi_norm_5)

    if use_entry_6:
        pprime_plot_6.set_data(psi_norm_6, pprime_6)
        ffprime_plot_6.set_data(psi_norm_6, ffprime_6)
        pressure_plot_6.set_data(psi_norm_6, pressure_6)
        psi_norm_plot_6.set_data(rho_tor_norm_6, psi_norm_6)

    if use_entry_7:
        pprime_plot_7.set_data(psi_norm_7, pprime_7)
        ffprime_plot_7.set_data(psi_norm_7, ffprime_7)
        pressure_plot_7.set_data(psi_norm_7, pressure_7)
        psi_norm_plot_7.set_data(rho_tor_norm_7, psi_norm_7)

    if use_entry_8:
        pprime_plot_8.set_data(psi_norm_8, pprime_8)
        ffprime_plot_8.set_data(psi_norm_8, ffprime_8)
        pressure_plot_8.set_data(psi_norm_8, pressure_8)
        psi_norm_plot_8.set_data(rho_tor_norm_8, psi_norm_8)

    if use_entry_9:
        pprime_plot_9.set_data(psi_norm_9, pprime_9)
        ffprime_plot_9.set_data(psi_norm_9, ffprime_9)
        pressure_plot_9.set_data(psi_norm_9, pressure_9)
        psi_norm_plot_9.set_data(rho_tor_norm_9, psi_norm_9)

    if use_entry_10:
        pprime_plot_10.set_data(psi_norm_10, pprime_10)
        ffprime_plot_10.set_data(psi_norm_10, ffprime_10)
        pressure_plot_10.set_data(psi_norm_10, pressure_10)
        psi_norm_plot_10.set_data(rho_tor_norm_10, psi_norm_10)

    

    # Update title with time
    ax1.set_title(f"$\psi_{{norm}}$ profile")    
    ax2.set_title(f"$j_{{\phi}}$ profile")    
    ax3.set_title(f"$\phi$ profile")
    ax4.set_title(f"q profiles")
    ax5.set_title(f"dpsi_drho_tor profile")
    ax6.set_title(f"gm1 profile")
    ax7.set_title(f"gm2 profile")
    ax8.set_title(f"gm3 profile")
    ax9.set_title(f"gm7 profile")
    ax10.set_title(f"dvolume_dpsi profile")

    # Set the main title
    fig.suptitle(f"1D profiles post inverse & post torax at time = {time[time_idx]:.2f}s", fontsize=14)

    # Formatting
    #ax1.set_aspect('equal')
    ax1.set_xlabel(r"$\rho_{tor,norm}$")
    ax1.set_ylabel(r"$\psi_{norm}$")
    #ax1.set_xlim(0, 1)
    #ax1.set_ylim(-1, 1)
    ax1.relim()
    ax1.autoscale_view(True, True, True)
    ax1.legend(fontsize=7, loc="best", frameon=True)
    

    #ax2.set_aspect('equal')
    ax2.set_xlabel(r"$\rho_{tor,norm}$")
    ax2.set_ylabel(r"$j_{{\phi}}$ (MA m$^{-2}$)")
    #ax2.set_xlim(0, 1)
    #ax2.set_ylim(-1, 1)
    ax2.relim()
    ax2.autoscale_view(True, True, True)
    ax2.legend(fontsize=7, loc="best", frameon=True)


    #ax3.set_aspect('equal')
    ax3.set_xlabel(r"$\rho_{tor,norm}$")
    ax3.set_ylabel(r"$\phi$ (Wb)")
    #ax3.set_xlim(0, 1)
    #ax3.set_ylim(-1, 1)
    ax3.relim()
    ax3.autoscale_view(True, True, True)
    ax3.legend(fontsize=7, loc="best", frameon=True)


    #ax4.set_aspect('equal')
    ax4.set_xlabel(r"$\rho_{tor,norm}$")
    ax4.set_ylabel(r"q")
    #ax4.set_xlim(0, 1)
    #ax4.set_ylim(-1, 1)
    ax4.relim()
    ax4.autoscale_view(True, True, True)
    ax4.legend(fontsize=7, loc="best", frameon=True)


    #ax5.set_aspect('equal')
    ax5.set_xlabel(r"$\rho_{tor,norm}$")
    ax5.set_ylabel(r"$\frac{d\psi}{d\rho_{tor}}$")
    #ax5.set_xlim(0, 1)
    #ax5.set_ylim(-1, 1)
    ax5.relim()
    ax5.autoscale_view(True, True, True)
    ax5.legend(fontsize=7, loc="best", frameon=True)

    ax6.set_xlabel(r"$\rho_{tor,norm}$")
    ax6.set_ylabel(r"gm1")
    ax6.relim()
    ax6.autoscale_view(True, True, True)
    ax6.legend(fontsize=7, loc="best", frameon=True)

    ax7.set_xlabel(r"$\rho_{tor,norm}$")
    ax7.set_ylabel(r"gm2")
    ax7.relim()
    ax7.autoscale_view(True, True, True)
    ax7.legend(fontsize=7, loc="best", frameon=True)

    ax8.set_xlabel(r"$\rho_{tor,norm}$")
    ax8.set_ylabel(r"gm3")
    ax8.relim()
    ax8.autoscale_view(True, True, True)
    ax8.legend(fontsize=7, loc="best", frameon=True)


    ax9.set_xlabel(r"$\rho_{tor,norm}$")
    ax9.set_ylabel(r"gm7")
    ax9.relim()
    ax9.autoscale_view(True, True, True)
    ax9.legend(fontsize=7, loc="best", frameon=True)


    ax10.set_xlabel(r"$\rho_{tor,norm}$")
    ax10.set_ylabel(r"$\frac{dV}{d\psi}$" )
    ax10.relim()
    ax10.autoscale_view(True, True, True)
    ax10.legend(fontsize=7, loc="best", frameon=True)


    
    # Refresh plot
    fig.canvas.draw_idle()

# Connect slider to update function
time_slider.on_changed(update_plot)

# Initial plot update
update_plot(0)

fig.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

if use_entry_0:
    entry_0.close()
if use_entry_1:
    entry_1.close()
if use_entry_2:
    entry_2.close()
if use_entry_3:
    entry_3.close()
if use_entry_4:
    entry_4.close()
if use_entry_5:
    entry_5.close()
if use_entry_6:
    entry_7.close()
if use_entry_7:
    entry_7.close()
if use_entry_8:
    entry_8.close()
if use_entry_9:
    entry_9.close()
if use_entry_10:
    entry_10.close()