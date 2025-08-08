import imas
import imaspy


# the purpose of this script is to combine IDS from two data entries to a new data entry

test_case = "test_23"
iteration = "iteration_10"

# Define data entry from input runs

entry_in_1 = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/"+test_case+"/initial_data","r")

entry_in_2 = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/"+test_case+"/"+iteration+"/post_torax","r")



# Read required IDSs from input data entries 


iron_core = entry_in_1.get("iron_core")
wall = entry_in_1.get("wall")
pf_passive = entry_in_1.get("pf_passive")
pf_active = entry_in_1.get("pf_active")
eq = entry_in_2.get("equilibrium")


# Create new data entry shot

entry_out = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/"+test_case+"/"+iteration+"/post_torax","w")


# Put the IDSs in the newly created data entry 

entry_out.put(iron_core)
entry_out.put(wall)
entry_out.put(pf_passive)
entry_out.put(pf_active)
entry_out.put(eq)


# Close all entries

entry_in_1.close()
entry_in_2.close()
entry_out.close()