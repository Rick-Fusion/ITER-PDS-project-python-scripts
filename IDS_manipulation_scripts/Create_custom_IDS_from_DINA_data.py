import imas
import numpy as np

""" The aim of this script is to create a custom data entry with IDSs from an existing data entry """


# define the time slices for which you want to extract data
time_slices_array = np.linspace(225,226,26)


# Open entry to read from
entry_in = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/3","r")
entry_out = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/6","w")

# Get and put the time slices that have no time dependence
wall = entry_in.get("wall")
iron_core = entry_in.get("iron_core")

entry_out.put(wall)
entry_out.put(iron_core)

for t in time_slices_array:
    eq = entry_in.get_slice("equilibrium", t, imas.ids_defs.CLOSEST_INTERP)
    pf_active = entry_in.get_slice("pf_active", t, imas.ids_defs.CLOSEST_INTERP)
    pf_passive = entry_in.get_slice("pf_passive", t, imas.ids_defs.CLOSEST_INTERP)

    entry_out.put_slice(eq)
    entry_out.put_slice(pf_active)
    entry_out.put_slice(pf_passive)

entry_in.close()
entry_out.close()