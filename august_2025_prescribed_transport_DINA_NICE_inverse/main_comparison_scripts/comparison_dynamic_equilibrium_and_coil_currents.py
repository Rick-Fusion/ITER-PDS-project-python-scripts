import imas
import imas.util
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patches as patches
import sys
from matplotlib.widgets import Slider



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



# This python script is dedicated to comparing an existing dynamic magnetic equilibrium of a tokamak plasma and the corresponding coil currents at a requested times by slider 
# with a newly computed equilibrium + coil currents with an inverse simulation.
# This script can for example be used to compare DINA results with NICE inverse results. This is useful for PDS deliverable 2 for the case of prescribed transport. 
#  - contours of the poloidal magnetic flux (psi) in a 2D plane at requested time with slider (for input data and output inverse data)
#  - poloidal magnetic flux map at requested time with slider (for input data and output inverse data)
#  - PF and CS coil cross-sections
#  - Vacuum vessel outline
#  - limiter wall + divertor outline
#  - PF and CS coil currents at the requested time


# For the plotting titles specify where input and output are coming from (Most likey input from DINA and output from NICE inverse)
input_data_source = "DINA"
output_data_source = "NICE inverse"



# Specify if the input and output magnetic equilibrium data is stored under profiles_2d (DINA) or ggd (NICE inverse)
equilibrium_in_data_stored = 'profiles_2d'      # either 'profiles_2d'  (DINA results) or  'ggd'  (NICE inverse results)

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



# SETTINGS FOR LOADING DATA ENTRY FOR PRESCRIBED TRANSPORT: DINA VS NICE INVERSE

SHOT = 105084
# SHOT = 105092

#  NOTE THAT FOR BOTH CASES A AND B THE DINA DATA IS PLOTTED WITH THE ORIGINAL VS COIL DATA. THIS MAKES SURE IT DOES NOT GIVE A DISTORTED VIEW ON THE OTHER COIL CURRENTS, SINCE THEY TAKE OVER WHAT VS IS NOT PROVIDING

# CASE = "A"          # NO REFERENCE | VS FIXED ON ZERO     (NO EXTERNAL DATA)
# CASE = "B"        # WITH REFERENCE | VS FIXED ON DINA   (EXTERNAL DATA)
# CASE = "C"          # WITH REFERENCE | VS FIXED ON ZERO  
CASE = "D"          # WITH REFERENCE | VS free ON DINA  


# time_slices_simulated = "interesting"
time_slices_simulated = "evenly_spaced"




### SHOT 105084 ###
if SHOT == 105084:


    if CASE == "A":

        if time_slices_simulated == "interesting":
        # CASE A: interesting time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - No reference currents - VS coil currents fixed to zero"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_A_interesting_time_slices", "r")

        elif time_slices_simulated == "evenly_spaced":
        # CASE A: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - No reference currents - VS coil currents fixed to zero"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_A_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()

    
    elif CASE == "B":
        
        if time_slices_simulated == "interesting":
        # CASE B: interesting time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents fixed to DINA data"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_B_interesting_time_slices", "r")

        elif time_slices_simulated == "evenly_spaced":
        # CASE B: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents fixed to DINA data"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_B_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()
    

    elif CASE == "C":
        
        # if time_slices_simulated == "interesting":
        # CASE C: interesting time slices

        if time_slices_simulated == "evenly_spaced":
        # CASE C: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents fixed to zero"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_C_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()


    elif CASE == "D":
        
        # if time_slices_simulated == "interesting":
        # CASE D: interesting time slices

        if time_slices_simulated == "evenly_spaced":
        # CASE D: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents free (DINA reference)"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_D_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()


    else: 
        print("The defined CASE is not a valid option")
        sys.exit()





### SHOT 105092 ###
elif SHOT == 105092:


    if CASE == "A":

        if time_slices_simulated == "interesting":
        # CASE A: interesting time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - No reference currents - VS coil currents fixed to zero"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_A_interesting_time_slices", "r")

        elif time_slices_simulated == "evenly_spaced":
        # CASE A: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - No reference currents - VS coil currents fixed to zero"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_A_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()

    
    elif CASE == "B":
        
        if time_slices_simulated == "interesting":
        # CASE B: interesting time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents fixed to DINA data"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_B_interesting_time_slices", "r")

        elif time_slices_simulated == "evenly_spaced":
        # CASE B: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents fixed to DINA data"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_B_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()
    
    
    elif CASE == "C":
        
        # if time_slices_simulated == "interesting":
        # CASE C: interesting time slices

        if time_slices_simulated == "evenly_spaced":
        # CASE C: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents fixed to zero"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_C_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()


    elif CASE == "D":
        
        # if time_slices_simulated == "interesting":
        # CASE D: interesting time slices

        if time_slices_simulated == "evenly_spaced":
        # CASE D: evenly spaced time slices
            global_title = rf"Prescribed transport: NICE inverse & DINA comparison - scenario: {SHOT} - With reference currents - VS coil currents free (DINA reference)"
            entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents", "r")
            entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/output_data_entries/NICE_inverse_prescribed_transport_cases/case_D_evenly_spaced_time_slices", "r")
        
        else:
            print("The defined time_slices_simulated is not a valid option")
            sys.exit()


    else: 
        print("The defined CASE is not a valid option")
        sys.exit()





else: 
    print("The defined SHOT is not a valid option")
    sys.exit()





entry_machine_description = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/3","r") # TODO: take geometry information from a general source? Machine description database?


# Inverse data can contain now multple time slices to plot with slider
# All the specific time slices that NICE inverse has simulated on will be filtered out from DINA data for direct comparison

eq_out = entry_out.get("equilibrium", lazy=False)
pf_active_out = entry_out.get("pf_active", lazy=False)
time_out_array = eq_out.time
number_of_time_slices = len(time_out_array)


print(f"The data that is loaded contains {number_of_time_slices} time slices.")


entry_in_sliced = imas.DBEntry("imas:memory?path=/","w")

for time in time_out_array:
    eq_in_slice = entry_in.get_slice("equilibrium", time, imas.ids_defs.CLOSEST_INTERP)
    pf_active_in_slice = entry_in.get_slice("pf_active", time, imas.ids_defs.CLOSEST_INTERP)
    entry_in_sliced.put_slice(eq_in_slice)
    entry_in_sliced.put_slice(pf_active_in_slice)


eq_in_sliced = entry_in_sliced.get("equilibrium", lazy=False)
pf_active_in_sliced = entry_in_sliced.get("pf_active", lazy=False)

wall = entry_machine_description.get("wall", lazy=True) 


### Assign required input and output data (equilibrium data and coil currents) ###


number_of_coils_in = pf_active_in_sliced.coil.size
number_of_coils_out = pf_active_out.coil.size

print(f"Amount of coils in input data = {number_of_coils_in}")
print(f"Amount of coils in output data = {number_of_coils_out}")



# Select which values of psi you want to plot a contour line (want to find the minimum and maximum for both equilibrium to decide ultimate max and min)
psi_min_in = float("inf")
psi_max_in = float("-inf")

psi_min_out = float("inf")
psi_max_out = float("-inf")

# When equilibrium data is under profiles_2d field: e.g. for DINA results
if equilibrium_in_data_stored == 'profiles_2d':
    for time_slice in eq_in_sliced.time_slice:
        psi_values = time_slice.profiles_2d[0].psi
        psi_min_in = min(psi_min_in, psi_values.min())
        psi_max_in = max(psi_max_in, psi_values.max())


# when equilibrium data is under ggd field: e.g. for NICE inverse results
if equilibrium_in_data_stored == 'ggd':
    for time_slice in eq_in_sliced.time_slice:
        psi_values = time_slice.ggd[0].psi[0].values
        psi_min_in = min(psi_min_in, psi_values.min())
        psi_max_in = max(psi_max_in, psi_values.max())



# When equilibrium data is under profiles_2d field: e.g. for DINA results
if equilibrium_out_data_stored == 'profiles_2d':
    for time_slice in eq_out.time_slice:
        psi_values = time_slice.profiles_2d[0].psi
        psi_min_out = min(psi_min_out, psi_values.min())
        psi_max_out = max(psi_max_out, psi_values.max())


# when equilibrium data is under ggd field: e.g. for NICE inverse results
if equilibrium_out_data_stored == 'ggd':
    for time_slice in eq_out.time_slice:
        psi_values = time_slice.ggd[0].psi[0].values
        psi_min_out = min(psi_min_out, psi_values.min())
        psi_max_out = max(psi_max_out, psi_values.max())




psi_min_global = min(psi_min_in, psi_min_out)
psi_max_global = max(psi_max_in, psi_max_out)


number_of_contours = 100                                    # Here specify the amount of psi iso contours that are visualised in the plot
contour_values_psi = np.linspace(psi_min_global, psi_max_global, number_of_contours) 


# input and output coil current data (put in np array)

coil_currents_in = np.zeros((number_of_coils_in, number_of_time_slices))
coil_currents_out = np.zeros((number_of_coils_out, number_of_time_slices))


for i in range(number_of_coils_in):
    for j in range(number_of_time_slices):
        coil_currents_in[i, j] = pf_active_in_sliced.coil[i].current.data[j]

for i in range(number_of_coils_out):
    for j in range(number_of_time_slices):
        coil_currents_out[i, j] = pf_active_out.coil[i].current.data[j]









# Create plot 

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(12,10))



# Formatting the plots (static)

ax1.set_aspect('equal')  # Equal aspect ratio
ax1.set_xlabel("R (m)")
ax1.set_ylabel("Z (m)")
ax1.set_xlim(0, 13)  # Adjust limits based on your data
ax1.set_ylim(-10, 10)



ax2.set_aspect('equal')  # Equal aspect ratio
ax2.set_xlabel("R (m)")
ax2.set_ylabel("Z (m)")
ax2.set_xlim(0, 13)  # Adjust limits based on your data
ax2.set_ylim(-10, 10)




# Set up a list of coil names (for plotting on x-axis of ax3)

coil_names = ['CS3U', 'CS2U', 'CS1U', 'CS1L', 'CS2L', 'CS3L',
              'PF1', 'PF2', 'PF3', 'PF4', 'PF5', 'PF6',
              'VSU', 'VSL']



ax3.set_xlabel('Coil Name')
ax3.set_ylabel('Current (kA)')
ax3.set_xticks(range(len(coil_names)))  # Set ticks at indices of coil_names
ax3.set_xticks([x - 0.5 for x in range(len(coil_names))], minor=True)  # Add gridlines between ticks
ax3.grid(which='minor', axis='x', linestyle='--', alpha=0.7)  # Gridlines between labels
ax3.grid(which='major', axis='y', linestyle='--', alpha=0.7)  # Normal gridlines for y-axis







# Color map & scatter plot settings

norm = mcolors.Normalize(vmin=psi_min_global, vmax=psi_max_global)          # Define normalization for consistent color mapping
#norm = mcolors.Normalize(vmin=-100, vmax=100)                              # Custom color scale

opaqueness_of_psi_scatter = 0                                               # at value 0.01 regions with fine mesh are visisble


# Create a hidden color mapping object for the color bar
sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
sm.set_array([])  # Required for the color bar


# Create the color bar on a separate axis, adjust its position using `cax`
cbar_ax = fig.add_axes([0.330, 0.15, 0.015, 0.7])  # [left, bottom, width, height]
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label("Poloidal magnetic flux (Wb)")







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


add_coil_rectangles(ax1, pf_active_in_sliced, number_of_coils_in)
add_coil_rectangles(ax2, pf_active_out, number_of_coils_out)






# The last two coils, Upper and Lower vertical Stabilization (VSL and VSU)coils, do not have rectangular shape, but the outline must be taken

VSU_outline_r = pf_active_in_sliced.coil[12].element[0].geometry.outline.r
VSU_outline_z = pf_active_in_sliced.coil[12].element[0].geometry.outline.z

VSL_outline_r = pf_active_in_sliced.coil[13].element[0].geometry.outline.r
VSL_outline_z = pf_active_in_sliced.coil[13].element[0].geometry.outline.z


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




### SETUP PLOTTING BASE ###




#  Initial plot data for initial time slice 
plasma_boundary_in, = ax1.plot([], [], '.', color="blue", alpha=0.8, markersize=2)
contour_lines_psi_in = None
scatter_psi_in = None


plasma_boundary_out, = ax2.plot([], [], '-', color="red", alpha=0.8, label='Inverse boundary') 
desired_plasma_boundary, = ax2.plot([], [], '.', color="blue", alpha=0.8, label='Desired boundary', markersize=2)
# desired_plasma_boundary, = ax2.plot([], [], '-', color="blue", alpha=0.5, label='Desired boundary')
contour_lines_psi_out = None
scatter_psi_out = None


coil_currents_in_plot = None
coil_currents_out_plot = None



# Display the poloidal magnetic flux (psi) at axis and boundary below the plot (Create text objects once)

psi_text_in = fig.text(
    0.5, -0.07,
    "",
    ha="center", va="top", transform=ax1.transAxes, fontsize=12, fontweight="bold"
)


psi_text_out = fig.text(
    0.5, -0.07,
    "",
    ha="center", va="top", transform=ax2.transAxes, fontsize=12, fontweight="bold"
)




# Add slider
ax_slider = plt.axes([0.2, -0.008, 0.6, 0.03], facecolor="lightgrey")
time_slider = Slider(ax_slider, "Time", 0, number_of_time_slices - 1, valinit=0, valstep=1)



# Handle key press events for slider control
fig.canvas.mpl_connect('key_press_event', on_key)


def update_plot(i):
    """Update the plot for a given time index i."""
    global contour_lines_psi_in, contour_lines_psi_out, coil_currents_in_plot, coil_currents_out_plot, scatter_psi_in, scatter_psi_out
    time_idx = int(i)



    ### QUANTITIES THAT CHANGE EVERY TIME SLICE ###




    # Interesting value of the poloidal magnetic flux (psi) are on the magnetic axis and plasma boundary - these will be displayed below the plots for the input and output data

    psi_axis_in = eq_in_sliced.time_slice[time_idx].global_quantities.psi_axis.value
    psi_boundary_in = eq_in_sliced.time_slice[time_idx].global_quantities.psi_boundary.value

    psi_axis_out = eq_out.time_slice[time_idx].global_quantities.psi_axis.value
    psi_boundary_out = eq_out.time_slice[time_idx].global_quantities.psi_boundary.value






    # Plasma boundary points of input equilibrium data - to plot the seperatrix in the input equilibrium contour plot

    plasma_boundary_r_in = eq_in_sliced.time_slice[time_idx].boundary.outline.r
    plasma_boundary_z_in = eq_in_sliced.time_slice[time_idx].boundary.outline.z






    # Plasma boundary points of output equilibrium data - to plot the seperatrix in the output equilibrium contour plot

    plasma_boundary_r_out = eq_out.time_slice[time_idx].boundary.outline.r
    plasma_boundary_z_out = eq_out.time_slice[time_idx].boundary.outline.z

 




    # Poloidal magnetic flux (psi) data from input equilibrium


    # When equilibrium data is under profiles_2d field: e.g. for DINA results

    if equilibrium_in_data_stored == 'profiles_2d':
        r_in = eq_in_sliced.time_slice[time_idx].profiles_2d[0].r
        z_in = eq_in_sliced.time_slice[time_idx].profiles_2d[0].z
        psi_in = eq_in_sliced.time_slice[time_idx].profiles_2d[0].psi

    # when equilibrium data is under ggd field: e.g. for NICE inverse results
    if equilibrium_in_data_stored == 'ggd':
        r_in = eq_in_sliced.time_slice[time_idx].ggd[0].r[0].values
        z_in = eq_in_sliced.time_slice[time_idx].ggd[0].z[0].values
        psi_in = eq_in_sliced.time_slice[time_idx].ggd[0].psi[0].values






    # Poloidal magnetic flux (psi) data from output equilibrium


    # When equilibrium data is under profiles_2d field: e.g. for DINA results

    if equilibrium_out_data_stored == 'profiles_2d':
        r_out = eq_out.time_slice[time_idx].profiles_2d[0].r
        z_out = eq_out.time_slice[time_idx].profiles_2d[0].z
        psi_out = eq_out.time_slice[time_idx].profiles_2d[0].psi

    # when equilibrium data is under ggd field: e.g. for NICE inverse results
    if equilibrium_out_data_stored == 'ggd':
        r_out = eq_out.time_slice[time_idx].ggd[0].r[0].values
        z_out = eq_out.time_slice[time_idx].ggd[0].z[0].values
        psi_out = eq_out.time_slice[time_idx].ggd[0].psi[0].values





    # Remove old values

    if contour_lines_psi_in is not None:
        for coll in contour_lines_psi_in.collections:
            coll.remove()

    if contour_lines_psi_out is not None:
        for coll in contour_lines_psi_out.collections:
            coll.remove()

    if coil_currents_in_plot is not None:
        coil_currents_in_plot.remove()

    if coil_currents_out_plot is not None:
        coil_currents_out_plot.remove()


    if scatter_psi_in is not None:
        scatter_psi_in.remove()

    if scatter_psi_out is not None:
        scatter_psi_out.remove()





    # Plot updated values


    plasma_boundary_in.set_data(plasma_boundary_r_in, plasma_boundary_z_in)
    contour_lines_psi_in = ax1.tricontour(r_in.flatten(), z_in.flatten(), psi_in.flatten(), levels = contour_values_psi, linewidths=1, cmap="viridis_r")
    scatter_psi_in = ax1.scatter(r_in, z_in, c=psi_in, cmap="viridis_r", alpha=opaqueness_of_psi_scatter, norm=norm)


    plasma_boundary_out.set_data(plasma_boundary_r_out, plasma_boundary_z_out)
    desired_plasma_boundary.set_data(plasma_boundary_r_in, plasma_boundary_z_in)
    contour_lines_psi_out = ax2.tricontour(r_out.flatten(), z_out.flatten(), psi_out.flatten(), levels = contour_values_psi, linewidths=1, cmap="viridis_r")
    scatter_psi_out = ax2.scatter(r_out, z_out, c=psi_out, cmap="viridis_r", alpha=opaqueness_of_psi_scatter, norm=norm)


    coil_currents_in_plot = ax3.scatter(coil_names, coil_currents_in[:,time_idx]/1e3, color='blue', s=20, label=f'Coil currents input (from {input_data_source})')       # From input pf_active IDS 
    coil_currents_out_plot = ax3.scatter(coil_names, coil_currents_out[:,time_idx]/1e3, color='orange', s=20, label=f'Coil currents output (from {output_data_source})')



    ax1.set_title(
        f"Input equilibrium poloidal magnetic flux from {input_data_source} at time: {time_out_array[time_idx]:.4f} s",
          fontsize=12, fontweight="bold", pad = 20
        )
    
    ax2.set_title(
        f"Output equilibrium poloidal magnetic flux from {output_data_source} at time: {time_out_array[time_idx]:.4f} s", 
        fontsize=12, fontweight="bold", pad = 20
        )
    
    ax3.set_title(
        f'in- and ouput coil currents at time: {time_out_array[time_idx]:.4f} s', fontsize=12, fontweight="bold"
        )


    psi_text_in.set_text(rf"Input:  $ψ_{{axis}}$={psi_axis_in:.3e} Wb,  $ψ_{{boundary}}$={psi_boundary_in:.4e} Wb")
    psi_text_out.set_text(rf"Output:  $ψ_{{axis}}$={psi_axis_out:.3e} Wb,  $ψ_{{boundary}}$={psi_boundary_out:.4e} Wb")
    


    # Refresh plot
    fig.canvas.draw_idle()




# Connect slider to update function
time_slider.on_changed(update_plot)

# Initial plot update
update_plot(0)


# Show legends
ax2.legend()
ax3.legend()

fig.suptitle(f"{global_title}", fontsize=16)


plt.tight_layout()
plt.show()