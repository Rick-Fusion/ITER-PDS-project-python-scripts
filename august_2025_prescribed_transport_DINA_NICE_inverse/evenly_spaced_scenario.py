import imas
import numpy as np
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.interpolate import interp1d as interp1
from imas import DBEntry, IDSFactory, convert_ids
from imas.ids_defs import CLOSEST_INTERP
import sys



N_TIMESLICES = 51               # amount of interesting time slices

SHOT = 105084
# SHOT = 105092



# Load the entry from which you want to take the time slices


if SHOT == 105084:

    # SHOT 105084
    #  entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_correct_boundary_and_profiles_zero_VS_coil_currents","r")
    entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_correct_boundary_and_profiles_original_VS_coil_currents","r")


    #  entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_zero_VS_coil_currents","w")
    entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents","w")



elif SHOT == 105092:
    # SHOT 105092
    # entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_correct_boundary_and_profiles_zero_VS_coil_currents","r")
    entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_correct_boundary_and_profiles_original_VS_coil_currents","r")

    # entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_zero_VS_coil_currents","w")
    entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_evenly_spaced_time_slices_correct_boundary_and_profiles_original_VS_coil_currents","w")



else: 
    print("The defined SHOT is not a valid option")
    sys.exit()



eq_load = entry_load.get("equilibrium", lazy=True)

time_array = eq_load.time


if SHOT == 105084:
    time_max = 276

elif SHOT == 105092:
    time_max = 166

time_array_evenly_spaced = np.linspace(time_array[0], time_max, N_TIMESLICES)

print(time_array_evenly_spaced)





wall = entry_load.get("wall", autoconvert=False)
iron_core = entry_load.get("iron_core", autoconvert=False)

skipped = []

for time in time_array_evenly_spaced:
    eq_load_slice = entry_load.get_slice("equilibrium", time, imas.ids_defs.CLOSEST_INTERP, autoconvert=False)
    pf_active_load_slice = entry_load.get_slice("pf_active", time, imas.ids_defs.CLOSEST_INTERP, autoconvert=False)
    pf_passive_load_slice = entry_load.get_slice("pf_passive", time, imas.ids_defs.CLOSEST_INTERP, autoconvert=False)

    # check if the slice is valid by looking if there is a boundary defined --> if not it will skip to the next slices
    if not eq_load_slice.time_slice[0].boundary.outline.has_value:
        skipped.append(time)
        continue


    entry_write.put_slice(eq_load_slice)
    entry_write.put_slice(pf_active_load_slice)
    entry_write.put_slice(pf_passive_load_slice)

entry_write.put(wall)
entry_write.put(iron_core)



entry_load.close()
entry_write.close()

print(rf"Writing slices to new DBEntry done! {len(skipped)} time slices were skipped, because they contained no boundary outline: {skipped}")