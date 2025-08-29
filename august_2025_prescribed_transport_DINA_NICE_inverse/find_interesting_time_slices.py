import imas
import numpy as np
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.interpolate import interp1d as interp1
from imas import DBEntry, IDSFactory, convert_ids
from imas.ids_defs import CLOSEST_INTERP



N_TIMESLICES = 51               # amount of interesting time slices



def find_interesting_time_slices(sm):
    t = sm.time
    # energy signal
    wth = sm.global_quantities.energy_thermal.value
    wbp = sm.global_quantities.energy_b_field_pol.value
    # for Bv
    ip = sm.global_quantities.ip.value
    li = sm.global_quantities.li.value
    betap = sm.global_quantities.beta_pol.value
    R = sm.boundary.magnetic_axis_r.value
    a = sm.boundary.minor_radius.value
    K = sm.boundary.elongation.value
    # indice_valid = (R>1) & (abs(ip) > 50e3)
    indice_valid = [i for i in range(len(R)) if R[i] > 1 and abs(ip[i]) > 50e3]
    R = max(np.array(R) + np.array([1]))
    # constante
    mu0 = 4 * np.pi * 1e-7
    # proxy for vertical magnetic field
    denom = 4 * np.pi * R * (8 * R / a / np.sqrt(K) + betap + li / 2 - 3 / 2)
    bv = mu0 * ip / denom
    # time derivative
    dwthdt = np.gradient(wth, t, edge_order=1)
    dwbpdt = np.gradient(wbp, t, edge_order=1)
    dbvdt = np.gradient(bv, t, edge_order=1)
    # control variable
    fwi = cumtrapz(t, abs(dwthdt) + abs(dwbpdt), initial=0)
    fwi = (fwi - min(fwi)) / (max(fwi) - min(fwi))
    fbvi = cumtrapz(t, abs(dbvdt), initial=0)
    fbvi = (fbvi - min(fbvi)) / (max(fbvi) - min(fbvi))
    # added to be strictely monotonic and have some points in flattop
    ft = (t - min(t)) / (max(t) - min(t))
    # f = (fwi + fbvi + ft) / 3
    f = ft
    # juste to take into account validity
    f = (f - min(f[indice_valid])) / (max(f[indice_valid] - min(f[indice_valid])))
    # time selection
    f_nearest = interp1(
        f[indice_valid],
        list(range(len(f[indice_valid]))),
        kind="linear",
    )
    indice_selected = sorted(
        list(set([int(idx) for idx in f_nearest(np.linspace(0, 1, N_TIMESLICES))]))
    )
    return indice_selected




# DINA DATA FROM WHICH THE SUMMARY IDS IS TAKEN TO DETERMINE INTERESTING TIME SLICES

# SHOT 105084
# entry_original_DINA = imas.DBEntry("imas:hdf5?path=/work/imas/shared/imasdb/ITER/3/105084/1","r")


# SHOT 105092
entry_original_DINA = imas.DBEntry("imas:mdsplus?path=/home/ITER/dubrovm/public/imasdb/iter/3/105092/1","r")


summary = entry_original_DINA.get("summary", autoconvert=False)
time_array = summary.time.value

interesting_time_slices = find_interesting_time_slices(summary)

print(interesting_time_slices)
print(time_array[interesting_time_slices])




# LOAD the entry from which you want to take the interesting time slices

# SHOT 105084
# entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_correct_boundary_and_profiles_zero_VS_coil_currents","r")
# entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_correct_boundary_and_profiles_original_VS_coil_currents","r")


# SHOT 105092
# entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_correct_boundary_and_profiles_zero_VS_coil_currents","r")
entry_load = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_correct_boundary_and_profiles_original_VS_coil_currents","r")






# WRITE to the entry where you want to store the interesting time slices

# SHOT 105084
# entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_zero_VS_coil_currents","w")
# entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/input_data_entries/DINA_data_3_105084_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_original_VS_coil_currents","w")


# SHOT 105092
# entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_zero_VS_coil_currents","w")
entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/input_data_entries/DINA_data_3_105092_1_DDv4.0.0_interesting_time_slices_correct_boundary_and_profiles_original_VS_coil_currents","w")






wall = entry_load.get("wall", autoconvert=False)
iron_core = entry_load.get("iron_core", autoconvert=False)

skipped = []

for time in time_array[interesting_time_slices]:
    eq_load_interesting_slice = entry_load.get_slice("equilibrium", time, imas.ids_defs.CLOSEST_INTERP, autoconvert=False)
    pf_active_load_interesting_slice = entry_load.get_slice("pf_active", time, imas.ids_defs.CLOSEST_INTERP, autoconvert=False)
    pf_passive_load_interesting_slice = entry_load.get_slice("pf_passive", time, imas.ids_defs.CLOSEST_INTERP, autoconvert=False)

    # check if the slice is valid by looking if there is a boundary defined --> if not it will skip to the next slices
    if not eq_load_interesting_slice.time_slice[0].boundary.outline.has_value:
        skipped.append(time)
        continue


    entry_write.put_slice(eq_load_interesting_slice)
    entry_write.put_slice(pf_active_load_interesting_slice)
    entry_write.put_slice(pf_passive_load_interesting_slice)

entry_write.put(wall)
entry_write.put(iron_core)



entry_original_DINA.close()
entry_load.close()
entry_write.close()

print(rf"Writing interesting slices to new DBEntry done! {len(skipped)} time slices were skipped, because they contained no boundary outline: {skipped}")