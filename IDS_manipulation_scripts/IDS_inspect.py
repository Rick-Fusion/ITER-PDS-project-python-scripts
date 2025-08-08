import imaspy
import imaspy.util
import matplotlib
import matplotlib.pyplot as plt


entry = imaspy.DBEntry(imaspy.ids_defs.HDF5_BACKEND, "ITER", 666666, 1, "vanschr", data_version="4")	#  define input data entry

entry.open()            # open input data entry

eq = entry.get("equilibrium")

imaspy.util.inspect(eq, hide_empty_nodes=True)		# inspect function


#imaspy.util.print_tree(eq....)				# print_tree function (only to be used for substructures of IDS)
