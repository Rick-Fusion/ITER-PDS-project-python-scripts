import imas
import numpy as np



entry_read_1 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/DINA_data_3_105084_1_DDv4.0.0_correct_boundary_correct_profiles","r")
# entry_read_2 = imas.DBEntry("imas:hdf5?path=/work/imas/shared/imasdb/ITER/3/105084/1","r")
entry_read_2 = imas.DBEntry("imas:mdsplus?path=/home/ITER/dubrovm/public/imasdb/iter/3/105092/1","r")


pf_active_1 = entry_read_1.get("pf_active")
eq = entry_read_1.get("equilibrium")
wall = entry_read_1.get("wall")
iron_core = entry_read_1.get("iron_core")
pf_passive = entry_read_1.get("pf_passive")         # pf_passive from DINA data is not supported by NICE inverse

pf_active_2 = entry_read_2.get("pf_active")




# for coil_idx in [12, 13]:
#     pf_active_1.coil[coil_idx].element[0] = pf_active_2.coil[coil_idx].element[0]


# # Set VS coil currents to zero
# pf_active_1.coil[12].current.data.value = np.zeros(len(pf_active_1.time))
# pf_active_1.coil[13].current.data.value = np.zeros(len(pf_active_1.time))


# Replace the value back to the original DINA data value
pf_active_1.coil[12].current.data.value = pf_active_2.coil[12].current.data.value
pf_active_1.coil[13].current.data.value = pf_active_2.coil[13].current.data.value


entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/DINA_data_3_105084_1_DDv4.0.0_correct_boundary_and_profiles_original_VS_coil_current","w")



entry_write.put(pf_active_1)
entry_write.put(eq)
entry_write.put(wall)
entry_write.put(iron_core)
entry_write.put(pf_passive)

entry_read_1.close()
# entry_read_2.close()
entry_write.close()