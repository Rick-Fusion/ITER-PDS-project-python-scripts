import imas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.ticker import FormatStrFormatter
matplotlib.use("TkAgg")

""" The aim of this plotting script is to visualise 0D profiles in general over time for multiple iterations resulting from the torax ouput or from any point in the NICE inverse + TORAX coupling scheme """

""" 
The 0D profiles that are currently plotted in this script include: 
- Total plasmma current ip


Currently there are two cases that occur for running this script and each time I have to adapt it to make it work. I will make an effort to make it more user friendly
for myself and perhaps others in the future. First I will identify which cases exist now:

- case 1: The initial data (entry_0) contains much more timeslices then the consecutive data entries. This is the case for test_1 where the original DINA data is used from vanschr 4/666666/3.
    - have two options: either load all the data and scale the plot for x and y axis (could be annoying) or selectively load only relevant time slices from entry_0 to compare
    - I find the latter nicer, since it now automatically scales the plot and you know the data is in the right time interval

- case 2: The initial data is already selectively extracted from the original data. This is done to have a custom time slice spacing and time interval for the equilibrium data
    - In this case the time interval is very similar for all data, so the data can be loaded as normal.

    
For case 1 you could choose to load all data for or selectively load data and for case 2 you would like to load all data. 

The choice is represented in the variable: selectively_load_data
- if True --> the initial_data is selected by taking time slices bases on the time_1 array that contains the time information from the data defined in entry_1
- if False --> all the data is loaded from the intiial_data and plotting must be scaled appropriately

"""

# Choice to selectively load data from initial_data based on entry_1 time (This contains the time interval based on simulation time)
selectively_load_data = False



#open the Data Entry
entry_0 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/initial_data", "r")
entry_1 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_1/post_torax", "r")
entry_2 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_2/post_torax", "r")
entry_3 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_3/post_torax", "r")
entry_4 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_4/post_torax", "r")
entry_5 = imas.DBEntry("imas:hdf5?path=/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme/test_8/iteration_5/post_torax", "r")


# depending on choice of fully loading intiial data or selectively, here the data is loaded fully immediately since it does not need time_1
if not selectively_load_data:
    eq_0 = entry_0.get("equilibrium", lazy=True)
    time_0 = eq_0.time


# get the full equilibrium IDS for all the iteration for 1 and higher
eq_1 = entry_1.get("equilibrium", lazy=True)
eq_2 = entry_2.get("equilibrium", lazy=True)
eq_3 = entry_3.get("equilibrium", lazy=True)
eq_4 = entry_4.get("equilibrium", lazy=True)
eq_5 = entry_5.get("equilibrium", lazy=True)

# Time data
time_1 = eq_1.time
time_2 = eq_2.time
time_3 = eq_3.time
time_4 = eq_4.time
time_5 = eq_5.time




# use this if you want to selectively load the initial data
# For entry_0 (the initial data), the times could not match, since only a certain interval is taken
# Therefore create the eq_0 for which time times do match with the other data entries by using time_1


if selectively_load_data:
    entry_0_mem = imas.DBEntry("imas:memory?path=/", "w")

    for time in time_1:
        eq_0_slice = entry_0.get_slice("equilibrium", time, imas.ids_defs.CLOSEST_INTERP)
        entry_0_mem.put_slice(eq_0_slice)

    eq_0 = entry_0_mem.get("equilibrium", lazy=True)

    time_0 = eq_0.time



print(f"## Time slices in eq_0 = {len(time_0)}##")
print(f"## Time slices in eq_1 = {len(time_1)}##")
print(f"## Time slices in eq_2 = {len(time_2)}##")
print(f"## Time slices in eq_3 = {len(time_3)}##")
print(f"## Time slices in eq_4 = {len(time_4)}##")
print(f"## Time slices in eq_5 = {len(time_5)}##")


# For each iteration, create an array of Ip, since now this is hidden in each time slice in the IDS 
# initialise ip arrays for each iteration
ip_0 = np.zeros(len(time_0))
ip_1 = np.zeros(len(time_1))
ip_2 = np.zeros(len(time_2))
ip_3 = np.zeros(len(time_3))
ip_4 = np.zeros(len(time_4))
ip_5 = np.zeros(len(time_5))


for i in range(0, len(time_0)):
    ip_0[i] = eq_0.time_slice[i].global_quantities.ip
    
for i in range(0,len(time_1)):   
    ip_1[i] = eq_1.time_slice[i].global_quantities.ip

for i in range(0,len(time_2)): 
    ip_2[i] = eq_2.time_slice[i].global_quantities.ip

for i in range(0,len(time_3)): 
    ip_3[i] = eq_3.time_slice[i].global_quantities.ip

for i in range(0,len(time_4)): 
    ip_4[i] = eq_4.time_slice[i].global_quantities.ip

for i in range(0,len(time_5)): 
    ip_5[i] = eq_5.time_slice[i].global_quantities.ip



# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(12, 10))

ax1.set_title(rf"Plasma current ($I_p$) against time for multiple iterations")
ax1.plot(time_0, ip_0/1e6, '.', color="green", label="ip initial_data")
ax1.plot(time_1, ip_1/1e6, '.', color="red", label="ip it=1")
ax1.plot(time_2, ip_2/1e6, '.', color="blue", label="ip it=2")
ax1.plot(time_3, ip_3/1e6, '.', color="magenta", label="ip it=3")
ax1.plot(time_4, ip_4/1e6, '.', color="yellow", label="ip it=4")
ax1.plot(time_5, ip_5/1e6, '.', color="black", label="ip it=5")

print("Difference between initial and last 5th iteration is:", ip_5-ip_0)

# Formatting
#ax1.set_aspect('equal')
ax1.set_xlabel(r"time (s)")
ax1.set_ylabel(r"$I_p (MA)$")
#ax1.set_xlim(225, 225.1)
#ax1.set_ylim(-15.01, -14.99)
ax1.relim()
ax1.autoscale_view(True, True, True)
ax1.legend(fontsize=12, loc="best", frameon=True)

plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
plt.tight_layout()
plt.show()


#entry_0_mem.close()
entry_0.close()
entry_1.close()
entry_2.close()
entry_3.close()
#entry_4.close()
#entry_5.close()
