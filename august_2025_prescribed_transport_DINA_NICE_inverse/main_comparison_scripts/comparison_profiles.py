import imas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import math
from scipy.interpolate import interp1d
import sys


matplotlib.use("TkAgg")


# function for controlling slider with left and right keys

def on_key(event):
    """Handle key presses to adjust the slider value."""
    if event.key == 'right':
        # Increase slider value by 1 (move right)
        slider.set_val(min(slider.val + 1, slider.valmax))
    elif event.key == 'left':
        # Decrease slider value by 1 (move left)
        slider.set_val(max(slider.val - 1, slider.valmin))


def get_nested_attr(obj, path):
    """Traverse nested attributes given a dotted path string."""
    for part in path.split("."):
        obj = getattr(obj, part)
    return obj


### THIS SCRIPT IS DELIBERATELY MADE FOR DINA DATA COMPARISON WITH NICE INVERSE SIMULATIONS FOR THE ITER SCENARIOS FOR THE PDS PROJECT ###
### DUE TO SOME DIFFERENCES IN IDS STRUCTURE THE DATA MUST BE HANDLED AND COMPARED IN A UNIQUE MANNER ###
### IT IS TAKEN INTO ACCOUNT THAT THE PSI PROFILE FROM DINA IS +1 IN LEGNTH COMPARED TO OTHER PROFILES FROM THE DINA EQUILIBRIUM IDS. THIS IS DONE TO MATCH THE PSI AT THE BOUNDARY ###
### IF INPUT DINA DATA IS USED WHERE THIS IS NOT THE CASE, CARE SHOULD BE TAKEN IN USING THIS SCRIPT ###



# Loading IDSs from the desired location (input and output)



# SETTINGS FOR LOADING DATA ENTRY FOR PRESCRIBED TRANSPORT: DINA VS NICE INVERSE

SHOT = 105084
# SHOT = 105092

#  NOTE THAT FOR BOTH CASES A AND B THE DINA DATA IS PLOTTED WITH THE ORIGINAL VS COIL DATA. THIS MAKES SURE IT DOES NOT GIVE A DISTORTED VIEW ON THE OTHER COIL CURRENTS, SINCE THEY TAKE OVER WHAT VS IS NOT PROVIDING

CASE = "A"          # NO REFERENCE | VS FIXED ON ZERO     (NO EXTERNAL DATA)
# CASE = "B"        # WITH REFERENCE | VS FIXED ON DINA   (EXTERNAL DATA)
# CASE = "C"          # WITH REFERENCE | VS FIXED ON ZERO  
# CASE = "D"          # WITH REFERENCE | VS free ON DINA  


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





eq_inverse = entry_out.get("equilibrium", lazy=True, autoconvert=False)
pf_active_inverse = entry_out.get("pf_active", lazy=True)

time_array = eq_inverse.time
n_slices = len(time_array) 


print(f"The data that is loaded contains {n_slices} time slices.")


entry_in_sliced = imas.DBEntry("imas:memory?path=/","w")

for time in time_array:
    eq_dina_slice = entry_in.get_slice("equilibrium", time, imas.ids_defs.CLOSEST_INTERP)
    pf_active_dina_slice = entry_in.get_slice("pf_active", time, imas.ids_defs.CLOSEST_INTERP)
    entry_in_sliced.put_slice(eq_dina_slice)
    entry_in_sliced.put_slice(pf_active_dina_slice)


eq_dina_sliced = entry_in_sliced.get("equilibrium", lazy=False)
pf_active_dina_sliced = entry_in_sliced.get("pf_active", lazy=False)



# --- Define quantities ---
profiles_1d = {
    "pressure": ("pressure", "pressure", "p [Pa]"),
    "f": ("f", "f", "f [m T]"),
    "pprime": ("dpressure_dpsi", "dpressure_dpsi", "p' [Pa / Wb]"),
    "ffprime": ("f_df_dpsi", "f_df_dpsi", "ff' [ $(m T)^2$ / Wb]"),
    "poloidal mangetic flux": ("psi", "psi", "ψ [Wb]"),
    "current density": ("j_phi", "j_phi", "$j_{\phi}$ [A / $m^2$]"),
    "safety factor": ("q", "q", "q"),
    "toroidal magnetic flux": ("phi", "phi", "$\phi$ [Wb]"),
    "volume": ("volume", "volume", "V [$m^3$]"),
    "surface": ("surface", "surface", "surface [$m^2$]"),
    "area": ("area", "area", "area [$m^2$]"),
    
    
}

boundary_quantities = {
    "elongation (@ $ψ_{norm}$ ~ 0.995)": ("boundary.elongation", "profiles_1d.elongation", "$\kappa$"),   # special case (profile vs scalar)
    # "triangularity_lower (NICE) (@ $ψ_{norm}$ ~ 0.995)": ("boundary.triangularity", "profiles_1d.triangularity_lower", "$\delta_{lower}$"), 
    # "triangularity_upper (NICE) (@ $ψ_{norm}$ ~ 0.995)": ("boundary.triangularity", "profiles_1d.triangularity_upper", "$\delta_{upper}$"),
    "triangularity (DINA @ $ψ_{norm}$ ~ 0.995 vs NICE @ $ψ_{norm}$ = 1)": ("boundary.triangularity", "boundary.triangularity", "$\delta$"),
    "$ψ_{separatrix}$ (@ $ψ_{norm}$ = 1)": ("boundary.psi", "boundary.psi", "ψ [Wb]"),
}

global_quantities = {
    "beta_tor": ("beta_tor", "beta_tor", r"$\beta_{tor}$"),
    "beta_pol": ("beta_pol", "beta_pol", r"$\beta_{pol}$"),
    "internal inductance": ("li_3", "li_3", "li_3"),
    "Total plasma current": ("ip", "ip", "Ip [A]"),
    "magnetic_axis_r": ("magnetic_axis.r", "magnetic_axis.r", "r [m]"),
    "magnetic_axis_z": ("magnetic_axis.z", "magnetic_axis.z", "z [m]"),
}


# --- Helper: match psi lengths ---
def match_psi(psi_norm, values):
    if len(values) < len(psi_norm):
        return psi_norm[:len(values)], values
    return psi_norm, values



# Prepare figure layout
n_profile = len(profiles_1d)
n_boundary = len(boundary_quantities)
n_global = len(global_quantities)

n_cols = 5
n_rows = math.ceil((n_profile + n_boundary + n_global) / n_cols)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
axes = axes.flatten()






# --- Precompute boundary + global quantities (time-dependent) ---
boundary_data_dina = {k: [] for k in boundary_quantities}
boundary_data_inverse = {k: [] for k in boundary_quantities}
global_data_dina = {k: [] for k in global_quantities}
global_data_inverse = {k: [] for k in global_quantities}

for i, t in enumerate(time_array):
    eq_dina_slice = eq_dina_sliced.time_slice[i]
    eq_inv_slice = eq_inverse.time_slice[i]

    # Boundary quantities
    for name, (dina_field, inverse_field, y_axis_label) in boundary_quantities.items():
        # --- DINA: always scalar from boundary.* ---
        dina_val = get_nested_attr(eq_dina_slice, dina_field).value
        boundary_data_dina[name].append(dina_val)

        # --- NICE inverse ---
    
        if inverse_field.startswith("profiles_1d."):
            # profile quantity → interpolate at psi=0.995
            psi_inv = eq_inv_slice.profiles_1d.psi.value
            psi_inv_norm = (psi_inv - psi_inv[0]) / (psi_inv[-1] - psi_inv[0])
            vals = get_nested_attr(eq_inv_slice, inverse_field)

            # determine psi_norm value on DINA defined boundary (probably close to 0.995, but can also depend on time)   This is also only for x-point case where the boundary is not correct
            psi_dina = eq_dina_slice.profiles_1d.psi.value
            psi_dina_norm = (psi_dina - psi_dina[0]) / (psi_dina[-1] - psi_dina[0])
            psi_dina_norm_boundary = psi_dina_norm[-2]                  # TODO: COMPUTE PSI NORM OF ORIGINAL DINA BOUNDARY (99.5% FLUX SURFACE) SUCH THAT THE BOUNDARY QUANTITIES ARE ALSO COMPARED (but this is not available anymore when this is overwritten in DDv4 IDS and last psi profiles value is replaced instead of added --> PROBLEM!)

            print(rf"The value of psi_norm at the original boundary from DINA at time {t} s = {psi_dina_norm_boundary}")

            f = interp1d(psi_inv_norm, vals, bounds_error=False, fill_value="extrapolate")
            inv_val = f(psi_dina_norm_boundary)
        else:
            # scalar quantity (e.g. boundary.elongation)
            inv_val = get_nested_attr(eq_inv_slice, inverse_field).value

        boundary_data_inverse[name].append(inv_val)

    # Global quantities
    for name, (dina_field, inverse_field, y_axis_label) in global_quantities.items():
        dina_val = get_nested_attr(eq_dina_slice.global_quantities, dina_field).value
        inv_val = get_nested_attr(eq_inv_slice.global_quantities, inverse_field).value
        global_data_dina[name].append(dina_val)
        global_data_inverse[name].append(inv_val)


# --- Plot boundary + global (time-dependent) ---
plot_idx = 0

time_markers = {}  # store vertical line handles

for name, (dina_field, inverse_field, y_axis_label) in boundary_quantities.items():
    axes[plot_idx].plot(time_array, boundary_data_dina[name], "-o", markersize=2 , color="blue", label="DINA")
    axes[plot_idx].plot(time_array, boundary_data_inverse[name], "-o", markersize=2, color="orange", label="NICE")
    axes[plot_idx].set_title(name)
    axes[plot_idx].set_xlabel("time [s]")
    axes[plot_idx].set_ylabel(rf"{y_axis_label}")
    axes[plot_idx].legend()
    # add vertical line at first time index
    vline = axes[plot_idx].axvline(time_array[0], color="red", linestyle="--")
    time_markers[name] = vline
    plot_idx += 1

for name, (dina_field, inverse_field, y_axis_label) in global_quantities.items():
    axes[plot_idx].plot(time_array, global_data_dina[name], "-o", markersize=2, color="blue", label="DINA")
    axes[plot_idx].plot(time_array, global_data_inverse[name], "-o", markersize=2, color="orange", label="NICE")
    axes[plot_idx].set_title(name)
    axes[plot_idx].set_xlabel("time [s]")
    axes[plot_idx].set_ylabel(rf"{y_axis_label}")
    axes[plot_idx].legend()
    # add vertical line at first time index
    vline = axes[plot_idx].axvline(time_array[0], color="red", linestyle="--")
    time_markers[name] = vline
    plot_idx += 1


# --- Profiles with slider ---
profile_axes = axes[plot_idx:plot_idx+n_profile]

# Initial plot at time slice 0
lines = {}
for ax, (name, (dina_field, inverse_field, y_axis_label)) in zip(profile_axes, profiles_1d.items()):
    psi_dina = eq_dina_sliced.time_slice[0].profiles_1d.psi.value
    psi_dina_norm = (psi_dina - psi_dina[0]) / (psi_dina[-1] - psi_dina[0])
    dina_vals = getattr(eq_dina_sliced.time_slice[0].profiles_1d, dina_field).value
    psi_dina_plot, dina_vals = match_psi(psi_dina_norm, dina_vals)

    psi_inv = eq_inverse.time_slice[0].profiles_1d.psi.value
    psi_inv_norm = (psi_inv - psi_inv[0]) / (psi_inv[-1] - psi_inv[0])
    inverse_vals = getattr(eq_inverse.time_slice[0].profiles_1d, inverse_field).value
    psi_inv_plot, inverse_vals = match_psi(psi_inv_norm, inverse_vals)

    l1, = ax.plot(psi_dina_plot, dina_vals, "-o", linewidth = 0.8, markersize=1, color="blue", label="DINA")
    l2, = ax.plot(psi_inv_plot, inverse_vals, "-o", linewidth = 0.8, markersize=1, color="orange", label="NICE")
    ax.set_title(f"{name} (t={time_array[0]:.2f}s)")
    ax.set_xlabel("Normalised ψ")
    ax.set_ylabel(rf"{y_axis_label}")
    ax.legend()
    lines[name] = (ax, l1, l2)

# Slider
ax_slider = plt.axes([0.2, 0.01, 0.6, 0.03])
slider = Slider(ax_slider, 'Time index', 0, n_slices-1, valinit=0, valstep=1)


# Handle key press events for slider control
fig.canvas.mpl_connect('key_press_event', on_key)

def update(idx):
    idx = int(idx)
    t_val = time_array[idx]
    for name, (dina_field, inverse_field, y_axis_label) in profiles_1d.items():
        ax, l1, l2 = lines[name]

        psi_dina = eq_dina_sliced.time_slice[idx].profiles_1d.psi.value
        psi_dina_norm = (psi_dina - psi_dina[0]) / (psi_dina[-1] - psi_dina[0])
        dina_vals = getattr(eq_dina_sliced.time_slice[idx].profiles_1d, dina_field).value
        psi_dina_plot, dina_vals = match_psi(psi_dina_norm, dina_vals)

        psi_inv = eq_inverse.time_slice[idx].profiles_1d.psi.value
        psi_inv_norm = (psi_inv - psi_inv[0]) / (psi_inv[-1] - psi_inv[0])
        inverse_vals = getattr(eq_inverse.time_slice[idx].profiles_1d, inverse_field).value
        psi_inv_plot, inverse_vals = match_psi(psi_inv_norm, inverse_vals)

        l1.set_xdata(psi_dina_plot)
        l1.set_ydata(dina_vals)
        l2.set_xdata(psi_inv_plot)
        l2.set_ydata(inverse_vals)
        ax.set_title(f"{name} (t={time_array[idx]:.2f}s)")
        ax.relim()
        ax.autoscale_view()

    # --- Update vertical markers in time-dependent plots ---
    for name, vline in time_markers.items():
        vline.set_xdata([t_val, t_val])


    fig.canvas.draw_idle()

slider.on_changed(update)

fig.suptitle(f"{global_title}", fontsize=16)


plt.tight_layout(rect=[0,0.05,1,1], pad=3.0, w_pad=1.0, h_pad=4.0)
plt.show()