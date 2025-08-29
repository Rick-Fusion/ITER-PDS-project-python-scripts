import imas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import sys



matplotlib.use("TkAgg")


# Loading IDSs from the desired location (input and output)



# SETTINGS FOR LOADING DATA ENTRY FOR PRESCRIBED TRANSPORT: DINA VS NICE INVERSE

# SHOT = 105084
SHOT = 105092

#  NOTE THAT FOR BOTH CASES A AND B THE DINA DATA IS PLOTTED WITH THE ORIGINAL VS COIL DATA. THIS MAKES SURE IT DOES NOT GIVE A DISTORTED VIEW ON THE OTHER COIL CURRENTS, SINCE THEY TAKE OVER WHAT VS IS NOT PROVIDING

# CASE = "A"          # NO REFERENCE | VS FIXED ON ZERO     (NO EXTERNAL DATA)
# CASE = "B"        # WITH REFERENCE | VS FIXED ON DINA   (EXTERNAL DATA)
CASE = "C"          # WITH REFERENCE | VS FIXED ON ZERO  
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




pf_active_out = entry_out.get("pf_active", lazy=False)
time_out_array = pf_active_out.time
number_of_time_slices = len(time_out_array)


print(f"The data that is loaded contains {number_of_time_slices} time slices.")


entry_in_sliced = imas.DBEntry("imas:memory?path=/","w")

for time in time_out_array:
    pf_active_in_slice = entry_in.get_slice("pf_active", time, imas.ids_defs.CLOSEST_INTERP)
    entry_in_sliced.put_slice(pf_active_in_slice)


pf_active_in_sliced = entry_in_sliced.get("pf_active", lazy=False)



number_of_coils_in = pf_active_in_sliced.coil.size
number_of_coils_out = pf_active_out.coil.size

print(f"Amount of coils in input data = {number_of_coils_in}")
print(f"Amount of coils in output data = {number_of_coils_out}")



coil_names = ['CS3U', 'CS2U', 'CS1U', 'CS1L', 'CS2L', 'CS3L',
              'PF1', 'PF2', 'PF3', 'PF4', 'PF5', 'PF6',
              'VSU', 'VSL']


# input and output coil current data (put in np array)

coil_currents_in = np.zeros((number_of_coils_in, number_of_time_slices))
coil_currents_out = np.zeros((number_of_coils_out, number_of_time_slices))


for i in range(number_of_coils_in):
    for j in range(number_of_time_slices):
        coil_currents_in[i, j] = pf_active_in_sliced.coil[i].current.data[j]

for i in range(number_of_coils_out):
    for j in range(number_of_time_slices):
        coil_currents_out[i, j] = pf_active_out.coil[i].current.data[j]



n_rows = 4
n_cols = 4

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))

fig.suptitle(f"{global_title}", fontsize=16)

# Flatten axes to make it iterable
axes = axes.flatten()


for ax, name, i in zip(axes, coil_names, range(number_of_coils_in)):
    ax.plot(time_out_array, coil_currents_in[i, :]/1e3, "-o", markersize=2, color="blue", label="DINA")
    ax.plot(time_out_array, coil_currents_out[i, :]/1e3, "-o", markersize=2, color="orange", label="NICE")
    ax.set_title(name)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Current [kA]")
    ax.legend()
    # if name == 'VSU' or  name == 'VSL':                                           # useful if both input as output VS coil currents are near zero. (but input (e.g. from DINA) is mostly not zero, but input for NICE inverse is set to zero )
    #     ax.set_ylim(bottom=-0.1, top=0.1)  # force near-zero band
    
# Hide any unused axes if coil_names < n_rows*n_cols
for ax in axes[number_of_coils_in:]:
    ax.set_visible(False)

plt.tight_layout(pad=3.0, w_pad=1.0, h_pad=2.0)
plt.show()