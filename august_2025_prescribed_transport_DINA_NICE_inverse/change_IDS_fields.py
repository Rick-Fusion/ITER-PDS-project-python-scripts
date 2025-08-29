import imas
import numpy as np

# The purpose of this script is to be able to specify the coil currents in a pf_active IDS:
# - based on coil currents from another pf_active IDS (which cannot be used directly in NICE Inverse for example)
# - custom values

# this functionality can be exploited from the plasma shape editor, but this is for now a fast option for testing.


entry_read_1 = imas.DBEntry("imas:mdsplus?path=/home/ITER/dubrovm/public/imasdb/iter/3/105092/1","r")
# entry_read_1 = imas.DBEntry("imas:hdf5?path=/work/imas/shared/imasdb/ITER/3/105084/1","r")

entry_read_2 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/DINA_data_3_105092_1_DDv4.0.0","r")

entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/DINA_data_3_105092_1_DDv4.0.0_correct_boundary_and_profiles","w")


eq_read_1 = entry_read_1.get("equilibrium", autoconvert=False)

eq_read_2 = entry_read_2.get("equilibrium", autoconvert=False)
pf_active = entry_read_2.get("pf_active", autoconvert=False)
pf_passive = entry_read_2.get("pf_passive", autoconvert=False)
wall = entry_read_2.get("wall", autoconvert=False)
iron_core = entry_read_2.get("iron_core", autoconvert=False)



time_array = eq_read_2.time

for i in range(len(time_array)):
    eq_read_2.time_slice[i].boundary.outline.r.value = eq_read_1.time_slice[i].boundary_separatrix.outline.r.value
    eq_read_2.time_slice[i].boundary.outline.z.value = eq_read_1.time_slice[i].boundary_separatrix.outline.z.value
    eq_read_2.time_slice[i].boundary.psi.value = -1 * eq_read_1.time_slice[i].boundary_separatrix.psi.value

    eq_read_2.time_slice[i].global_quantities.psi_boundary.value = -1 * eq_read_1.time_slice[i].boundary_separatrix.psi.value

    # --- profiles_1d arrays ---
    psi_array = eq_read_2.time_slice[i].profiles_1d.psi.value
    pprime_array = eq_read_2.time_slice[i].profiles_1d.dpressure_dpsi.value
    ffprime_array = eq_read_2.time_slice[i].profiles_1d.f_df_dpsi.value
    pressure_array = eq_read_2.time_slice[i].profiles_1d.pressure.value

    psi_separatrix = -1 * eq_read_1.time_slice[i].boundary_separatrix.psi.value

    new_psi_array = np.append(psi_array, psi_separatrix)
    new_pprime_array = np.append(pprime_array, 0.0)
    new_ffprime_array = np.append(ffprime_array, 0.0)
    new_pressure_array = np.append(pressure_array, 0.0)

    # Assign back (IDSNumericArray will be rebuilt automatically)
    eq_read_2.time_slice[i].profiles_1d.psi = new_psi_array
    eq_read_2.time_slice[i].profiles_1d.dpressure_dpsi = new_pprime_array
    eq_read_2.time_slice[i].profiles_1d.f_df_dpsi = new_ffprime_array
    eq_read_2.time_slice[i].profiles_1d.pressure = new_pressure_array







entry_write.put(eq_read_2)
entry_write.put(pf_active)
entry_write.put(pf_passive)
entry_write.put(wall)
entry_write.put(iron_core)


entry_read_1.close()
entry_read_2.close()
entry_write.close()