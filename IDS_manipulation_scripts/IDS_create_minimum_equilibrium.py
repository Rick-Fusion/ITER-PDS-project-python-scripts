import imas
import imaspy

# The purpose of this script is to create a new DBEntry with copied over pf_active, pf_passive, wall and iron_core IDSs
# and create an empty equilibrium IDS and copy over only relevant fields for testing NICE inverse input requirements.



# Define DB entry to read data from and define new one
entry_in = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/301","r")
entry_out = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/302","w")



# Read IDSs from input entry that are used to set the output entry 
iron_core = entry_in.get("iron_core")
wall = entry_in.get("wall")
pf_passive = entry_in.get("pf_passive")
pf_active = entry_in.get("pf_active")

eq_in = entry_in.get("equilibrium")



# create empty equilibrium IDS  
eq_out = imaspy.IDSFactory().equilibrium()


# Set mandatory field for each IDS (is this really mandatory?)
eq_out.ids_properties.homogeneous_time = eq_in.ids_properties.homogeneous_time #Should be 0 or 1 ?
eq_out.ids_properties.comment = "custom filled empty IDS to test NICE inverse input requirements."



# Set relevant fields for NICE inverse input from input equilibrium IDS
eq_out.time = eq_in.time



#eq_out.time_slice = eq_in.time_slice
eq_out.time_slice.resize(1)

#eq_out.time_slice[0] = eq_in.time_slice[0]

eq_out.time_slice[0].boundary.outline = eq_in.time_slice[0].boundary.outline

eq_out.time_slice[0].global_quantities.ip = eq_in.time_slice[0].global_quantities.ip

eq_out.time_slice[0].profiles_1d.psi = eq_in.time_slice[0].profiles_1d.psi
eq_out.time_slice[0].profiles_1d.dpressure_dpsi = eq_in.time_slice[0].profiles_1d.dpressure_dpsi
eq_out.time_slice[0].profiles_1d.f_df_dpsi = eq_in.time_slice[0].profiles_1d.f_df_dpsi


eq_out.vacuum_toroidal_field = eq_in.vacuum_toroidal_field


# This is now the minimum necessary in the equilibrium IDS to be able to run NICE inverse !!




# Put the IDSs to the output data entry 

entry_out.put(iron_core)
entry_out.put(wall)
entry_out.put(pf_passive)
entry_out.put(pf_active)
entry_out.put(eq_out)


# close entries

entry_in.close()
entry_out.close()