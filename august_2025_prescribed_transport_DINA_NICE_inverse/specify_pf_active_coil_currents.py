import imas
import numpy as np

# The purpose of this script is to be able to specify the coil currents in a pf_active IDS:
# - based on coil currents from another pf_active IDS (which cannot be used directly in NICE Inverse for example)
# - custom values

# this functionality can be exploited from the plasma shape editor, but this is for now a fast option for testing.



entry_read_for_pf_active = imas.DBEntry("imas:hdf5?path=/work/imas/shared/imasdb/ITER/3/105084/1","r")
entry_read_DDv4 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/DINA_data_3_105084_1_DDv4.0.0_random_pf_active_wrong_boundary","r")


entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3.5MA_2.65T_H/DINA_data_3_105084_1_DDv4.0.0_DINA_pf_active_right_boundary","w")


pf_active_read = entry_read_for_pf_active.get("pf_active")

pf_active_random = entry_read_DDv4.get("pf_active")
eq = entry_read_DDv4.get("equilibrium")
wall = entry_read_DDv4.get("wall")
iron_core = entry_read_DDv4.get("iron_core")
pf_passive = entry_read_DDv4.get("pf_passive")



# Set the time of pf_active equal to time from which you read the pf_active coil currents
pf_active_random.time = pf_active_read.time



# Set the coil currents for all times equal to the ones from the desired data (not vor VS coils, since these will be set to zero)
for coil_index in range(14):
    pf_active_random.coil[coil_index].current.data.value = pf_active_read.coil[coil_index].current.data.value
    pf_active_random.coil[coil_index].voltage.data.value = pf_active_read.coil[coil_index].voltage.data.value
    pf_active_random.supply[coil_index].voltage.data.value = pf_active_read.supply[coil_index].voltage.data.value


# Set VS coil currents to zero
pf_active_random.coil[12].current.data.value = np.zeros(len(pf_active_read.time))
pf_active_random.coil[13].current.data.value = np.zeros(len(pf_active_read.time))



entry_write.put(pf_active_random)
entry_write.put(eq)
entry_write.put(wall)
entry_write.put(iron_core)
entry_write.put(pf_passive)

entry_read_DDv4.close()
entry_read_for_pf_active.close()
entry_write.close()