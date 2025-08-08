import imas
import imaspy
import imaspy.util


# The purpose of this script is to copy field from one IDS to another IDS


entry1 = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/3","r")
#entry2 = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/101","r")


eq1 = entry1.get("equilibrium", lazy=True)

time_array = eq1.time

wall = entry1.get("wall")
iron_core = entry1.get("iron_core")



with imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/4", "w") as db_entry:
    db_entry.put(wall)
    db_entry.put(iron_core)
    
    for time in time_array[0:10]:
        
        # Get all the time slices from old DB entry
        eq = entry1.get_slice("equilibrium", time, imaspy.ids_defs.CLOSEST_INTERP)    # Lazy loaded IDSs cannot be used in put or put_slice, and this more efficient than loading completely
        pf_active = entry1.get_slice("pf_active", time, imaspy.ids_defs.CLOSEST_INTERP)
        pf_passive = entry1.get_slice("pf_passive", time, imaspy.ids_defs.CLOSEST_INTERP)

        
        # Put all the time slices to new DB entry
        db_entry.put_slice(eq)
        db_entry.put_slice(pf_active)
        db_entry.put_slice(pf_passive)



entry1.close()
#entry2.close()