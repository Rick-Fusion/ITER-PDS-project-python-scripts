import imas
import numpy as np



entry_read_1 = imas.DBEntry("imas:mdsplus?path=/home/ITER/dubrovm/public/imasdb/iter/3/105092/1","r")
# entry_read_1 = imas.DBEntry("imas:hdf5?path=/work/imas/shared/imasdb/ITER/3/105084/1","r")
entry_read_2 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/3","r")

entry_write = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/pds_test_12_08_2025/pds/run/simulations/test_DINA_scenario_3MA_2.65T_H/DINA_data_3_105092_1_DDv3.39.0","w")


# eq = entry_read_1.get("equilibrium", autoconvert=False)
# pf_active = entry_read_1.get("pf_active", autoconvert=False)
# pf_passive = entry_read_1.get("pf_passive", autoconvert=False)
# wall = entry_read_1.get("wall", autoconvert=False)
# iron_core = entry_read_2.get("iron_core", autoconvert=False)

# !!! WHEN PUTTING IDSs MAKE SURE THAT THE IMAS MODULE LOADED IS IN THE VERSION YOU WANT THE IDSs TO BE LABELED (can cause issue if trying to perform correct DD conversions afterwards!!)

eq = entry_read_1.get("equilibrium")
pf_active = entry_read_1.get("pf_active")
pf_passive = entry_read_1.get("pf_passive")
wall = entry_read_1.get("wall")
iron_core = entry_read_2.get("iron_core")


entry_write.put(eq)
entry_write.put(pf_active)
entry_write.put(pf_passive)
entry_write.put(wall)
entry_write.put(iron_core)


entry_read_1.close()
entry_read_2.close()
entry_write.close()