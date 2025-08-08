import imas
import imaspy
import imaspy.util

# The purpose of this script is to compare two IDSs with eachother


entry1 = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/101","r")
entry2 = imaspy.DBEntry("imas:hdf5?path=/home/ITER/vanschr/public/imasdb/ITER/4/666666/201","r")


eq1 = entry1.get("equilibrium")
eq2 = entry2.get("equilibrium")



imaspy.util.inspect(eq1.vacuum_toroidal_field, hide_empty_nodes=True)
imaspy.util.inspect(eq2.vacuum_toroidal_field, hide_empty_nodes=True)

#entry1.close()
#entry2.close()