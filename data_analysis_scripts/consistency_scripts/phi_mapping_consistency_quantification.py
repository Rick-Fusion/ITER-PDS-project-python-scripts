import imas
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider
matplotlib.use("TkAgg")



""" The aim of this script is to be able to visualise data consistency of mapping in terms of RMSE with varying n_nice and n_torax in a 2D data map """



def load_equilibria_it1_post_inverse_torax(test_case_number, base_path):
    """
    Load IMAS equilibrium entries for a specific test case for iteration 1 post inverse and post torax data

    Args:
        test_case_number (int): Test case number (j).
        base_path (str): Base path where test folders are stored.

    Returns:
        List of equilibrium objects, with the initial data at index 0.
    """
    equilibria = []

    # Construct base test case path
    test_case_path = os.path.join(base_path, f"test_{test_case_number}")


    iteration_path_nice = os.path.join(test_case_path, f"iteration_1", "post_inverse")
    uri = f"imas:hdf5?path={iteration_path_nice}"
    entry = imas.DBEntry(uri, "r")
    eq = entry.get("equilibrium", lazy=True)
    equilibria.append(eq)
    
 
    iteration_path_torax = os.path.join(test_case_path, f"iteration_1", "post_torax")
    uri = f"imas:hdf5?path={iteration_path_torax}"
    entry = imas.DBEntry(uri, "r")
    eq = entry.get("equilibrium", lazy=True)
    equilibria.append(eq)

    return equilibria



def compute_rmse(array1, array2, x1=None, x2=None):
    """
    Computes RMSE between two 1D profiles, interpolating to their overlapping domain if needed.

    If the profiles and their grids match in shape and values, RMSE is directly computed.
    Otherwise, interpolation is used within the overlapping domain.

    Args:
        array1 (np.ndarray): First profile (e.g. phi_1).
        array2 (np.ndarray): Second profile (e.g. phi_2).
        x1 (np.ndarray, optional): Grid for array1 (e.g. rho_1).
        x2 (np.ndarray, optional): Grid for array2 (e.g. rho_2).

    Returns:
        float: RMSE within overlapping domain.

    Raises:
        ValueError: If no overlapping domain is found or if grid information is missing when needed.
    """
    array1 = np.asarray(array1)
    array2 = np.asarray(array2)

    # If arrays match in shape but x1/x2 are given, check if grids match
    if array1.shape == array2.shape:
        if x1 is not None and x2 is not None:
            x1 = np.asarray(x1)
            x2 = np.asarray(x2)

            if x1.shape != array1.shape or x2.shape != array2.shape:
                raise ValueError("x1 and array1 (or x2 and array2) must have the same shape.")

            if np.allclose(x1, x2):
                # Same grid: direct comparison
                return np.sqrt(np.mean((array1 - array2) ** 2))
            else:
                # Same shape but different grid: interpolate
                pass  # proceed to interpolation below
        else:
            # No x1/x2 provided, but same shape: direct comparison
            return np.sqrt(np.mean((array1 - array2) ** 2))

    # Interpolation case
    if x1 is None or x2 is None:
        raise ValueError("Shapes differ — must provide both x1 and x2 to enable interpolation and overlapping comparison.")

    x1 = np.asarray(x1)
    x2 = np.asarray(x2)

    if array1.shape != x1.shape:
        raise ValueError("array1 and x1 must have the same shape.")
    if array2.shape != x2.shape:
        raise ValueError("array2 and x2 must have the same shape.")

    # Determine overlapping range
    overlap_low = max(x1.min(), x2.min()) 
    overlap_high = min(x1.max(), x2.max())    

    if overlap_low >= overlap_high:
        raise ValueError("No overlapping domain for comparison.")

    # Create new grid within overlapping range
    num_new = min(len(x1), len(x2))
    new_grid = np.linspace(overlap_low, overlap_high, num_new)

    # Interpolate both to new grid
    interpolator1 = interp1d(x1, array1, kind='linear', bounds_error=True)
    interpolator2 = interp1d(x2, array2, kind='linear', bounds_error=True)

    array1_new = interpolator1(new_grid)
    array2_new = interpolator2(new_grid)

    # Calculate RMSE on overlapping domain
    rmse = np.sqrt(np.mean((array1_new - array2_new) ** 2))

    return rmse



# load data from all test cases


# Test case A
equilibria_n_NICE_50_n_TORAX_10 = load_equilibria_it1_post_inverse_torax(test_case_number = 10, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_50_n_TORAX_15 = load_equilibria_it1_post_inverse_torax(test_case_number = 9, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_50_n_TORAX_20 = load_equilibria_it1_post_inverse_torax(test_case_number = 8, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_50_n_TORAX_25 = load_equilibria_it1_post_inverse_torax(test_case_number = 7.2, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_50_n_TORAX_30 = load_equilibria_it1_post_inverse_torax(test_case_number = 12, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_50_n_TORAX_50 = load_equilibria_it1_post_inverse_torax(test_case_number = 13, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_50_n_TORAX_100 = load_equilibria_it1_post_inverse_torax(test_case_number = 11, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")


# Test case B (2 cases overlap with other two cases)
equilibria_n_NICE_25_n_TORAX_25 = load_equilibria_it1_post_inverse_torax(test_case_number = 17, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_100_n_TORAX_25 = load_equilibria_it1_post_inverse_torax(test_case_number = 14, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")


# Test case C
equilibria_n_NICE_1000_n_TORAX_10 = load_equilibria_it1_post_inverse_torax(test_case_number = 18, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_1000_n_TORAX_15 = load_equilibria_it1_post_inverse_torax(test_case_number = 23, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_1000_n_TORAX_20 = load_equilibria_it1_post_inverse_torax(test_case_number = 21, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_1000_n_TORAX_25 = load_equilibria_it1_post_inverse_torax(test_case_number = 15, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_1000_n_TORAX_30 = load_equilibria_it1_post_inverse_torax(test_case_number = 22, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_1000_n_TORAX_50 = load_equilibria_it1_post_inverse_torax(test_case_number = 20, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")
equilibria_n_NICE_1000_n_TORAX_100 = load_equilibria_it1_post_inverse_torax(test_case_number = 19, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")



# Load the data from equilibrium IDS
# Test case A
phi_list_n_NICE_50_n_TORAX_10 = [eq.time_slice[0].profiles_1d.phi for eq in equilibria_n_NICE_50_n_TORAX_10]
rho_tor_norm_list_n_NICE_50_n_TORAX_10 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_50_n_TORAX_10]

phi_list_n_NICE_50_n_TORAX_15 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_50_n_TORAX_15]
rho_tor_norm_list_n_NICE_50_n_TORAX_15 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_50_n_TORAX_15]

phi_list_n_NICE_50_n_TORAX_20 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_50_n_TORAX_20]
rho_tor_norm_list_n_NICE_50_n_TORAX_20 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_50_n_TORAX_20]

phi_list_n_NICE_50_n_TORAX_25 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_50_n_TORAX_25]
rho_tor_norm_list_n_NICE_50_n_TORAX_25 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_50_n_TORAX_25]

phi_list_n_NICE_50_n_TORAX_30 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_50_n_TORAX_30]
rho_tor_norm_list_n_NICE_50_n_TORAX_30 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_50_n_TORAX_30]

phi_list_n_NICE_50_n_TORAX_50 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_50_n_TORAX_50]
rho_tor_norm_list_n_NICE_50_n_TORAX_50 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_50_n_TORAX_50]

phi_list_n_NICE_50_n_TORAX_100 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_50_n_TORAX_100]
rho_tor_norm_list_n_NICE_50_n_TORAX_100 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_50_n_TORAX_100]


# Test case B
phi_list_n_NICE_25_n_TORAX_25 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_25_n_TORAX_25]
rho_tor_norm_list_n_NICE_25_n_TORAX_25 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_25_n_TORAX_25]

phi_list_n_NICE_100_n_TORAX_25 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_100_n_TORAX_25]
rho_tor_norm_list_n_NICE_100_n_TORAX_25 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_100_n_TORAX_25]


# Test case C
phi_list_n_NICE_1000_n_TORAX_10 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_1000_n_TORAX_10]
rho_tor_norm_list_n_NICE_1000_n_TORAX_10 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_1000_n_TORAX_10]

phi_list_n_NICE_1000_n_TORAX_15 = [eq.time_slice[0].profiles_1d.phi for eq in equilibria_n_NICE_1000_n_TORAX_15]
rho_tor_norm_list_n_NICE_1000_n_TORAX_15 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_1000_n_TORAX_15]

phi_list_n_NICE_1000_n_TORAX_20 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_1000_n_TORAX_20]
rho_tor_norm_list_n_NICE_1000_n_TORAX_20 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_1000_n_TORAX_20]

phi_list_n_NICE_1000_n_TORAX_25 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_1000_n_TORAX_25]
rho_tor_norm_list_n_NICE_1000_n_TORAX_25 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_1000_n_TORAX_25]

phi_list_n_NICE_1000_n_TORAX_30 = [eq.time_slice[0].profiles_1d.phi for eq in equilibria_n_NICE_1000_n_TORAX_30]
rho_tor_norm_list_n_NICE_1000_n_TORAX_30 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_1000_n_TORAX_30]

phi_list_n_NICE_1000_n_TORAX_50 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_1000_n_TORAX_50]
rho_tor_norm_list_n_NICE_1000_n_TORAX_50 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_1000_n_TORAX_50]

phi_list_n_NICE_1000_n_TORAX_100 = [eq.time_slice[0].profiles_1d.phi  for eq in equilibria_n_NICE_1000_n_TORAX_100]
rho_tor_norm_list_n_NICE_1000_n_TORAX_100 = [eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_n_NICE_1000_n_TORAX_100]



# determine rmse for each case now:

# Test case A
rmse_n_NICE_50_n_TORAX_10 = compute_rmse(phi_list_n_NICE_50_n_TORAX_10[0], phi_list_n_NICE_50_n_TORAX_10[1], rho_tor_norm_list_n_NICE_50_n_TORAX_10[0], rho_tor_norm_list_n_NICE_50_n_TORAX_10[1])

rmse_n_NICE_50_n_TORAX_15 = compute_rmse(phi_list_n_NICE_50_n_TORAX_15[0], phi_list_n_NICE_50_n_TORAX_15[1], rho_tor_norm_list_n_NICE_50_n_TORAX_15[0], rho_tor_norm_list_n_NICE_50_n_TORAX_15[1])

rmse_n_NICE_50_n_TORAX_20 = compute_rmse(phi_list_n_NICE_50_n_TORAX_20[0], phi_list_n_NICE_50_n_TORAX_20[1], rho_tor_norm_list_n_NICE_50_n_TORAX_20[0], rho_tor_norm_list_n_NICE_50_n_TORAX_20[1])

rmse_n_NICE_50_n_TORAX_25 = compute_rmse(phi_list_n_NICE_50_n_TORAX_25[0], phi_list_n_NICE_50_n_TORAX_25[1], rho_tor_norm_list_n_NICE_50_n_TORAX_25[0], rho_tor_norm_list_n_NICE_50_n_TORAX_25[1])

rmse_n_NICE_50_n_TORAX_30 = compute_rmse(phi_list_n_NICE_50_n_TORAX_30[0], phi_list_n_NICE_50_n_TORAX_30[1], rho_tor_norm_list_n_NICE_50_n_TORAX_30[0], rho_tor_norm_list_n_NICE_50_n_TORAX_30[1])

rmse_n_NICE_50_n_TORAX_50 = compute_rmse(phi_list_n_NICE_50_n_TORAX_50[0], phi_list_n_NICE_50_n_TORAX_50[1], rho_tor_norm_list_n_NICE_50_n_TORAX_50[0], rho_tor_norm_list_n_NICE_50_n_TORAX_50[1])

rmse_n_NICE_50_n_TORAX_100 = compute_rmse(phi_list_n_NICE_50_n_TORAX_100[0], phi_list_n_NICE_50_n_TORAX_100[1], rho_tor_norm_list_n_NICE_50_n_TORAX_100[0], rho_tor_norm_list_n_NICE_50_n_TORAX_100[1])


# Test case B
rmse_n_NICE_25_n_TORAX_25 = compute_rmse(phi_list_n_NICE_25_n_TORAX_25[0], phi_list_n_NICE_25_n_TORAX_25[1], rho_tor_norm_list_n_NICE_25_n_TORAX_25[0], rho_tor_norm_list_n_NICE_25_n_TORAX_25[1])

rmse_n_NICE_100_n_TORAX_25 = compute_rmse(phi_list_n_NICE_100_n_TORAX_25[0], phi_list_n_NICE_100_n_TORAX_25[1], rho_tor_norm_list_n_NICE_100_n_TORAX_25[0], rho_tor_norm_list_n_NICE_100_n_TORAX_25[1])


# Test case C
rmse_n_NICE_1000_n_TORAX_10 = compute_rmse(phi_list_n_NICE_1000_n_TORAX_10 [0], phi_list_n_NICE_1000_n_TORAX_10 [1], rho_tor_norm_list_n_NICE_1000_n_TORAX_10 [0], rho_tor_norm_list_n_NICE_1000_n_TORAX_10 [1])

rmse_n_NICE_1000_n_TORAX_15 = compute_rmse(phi_list_n_NICE_1000_n_TORAX_15[0], phi_list_n_NICE_1000_n_TORAX_15[1], rho_tor_norm_list_n_NICE_1000_n_TORAX_15[0], rho_tor_norm_list_n_NICE_1000_n_TORAX_15[1])

rmse_n_NICE_1000_n_TORAX_20 = compute_rmse(phi_list_n_NICE_1000_n_TORAX_20[0], phi_list_n_NICE_1000_n_TORAX_20[1], rho_tor_norm_list_n_NICE_1000_n_TORAX_20[0], rho_tor_norm_list_n_NICE_1000_n_TORAX_20[1])

rmse_n_NICE_1000_n_TORAX_25 = compute_rmse(phi_list_n_NICE_1000_n_TORAX_25[0], phi_list_n_NICE_1000_n_TORAX_25[1], rho_tor_norm_list_n_NICE_1000_n_TORAX_25[0], rho_tor_norm_list_n_NICE_1000_n_TORAX_25[1])

rmse_n_NICE_1000_n_TORAX_30 = compute_rmse(phi_list_n_NICE_1000_n_TORAX_30[0], phi_list_n_NICE_1000_n_TORAX_30[1], rho_tor_norm_list_n_NICE_1000_n_TORAX_30[0], rho_tor_norm_list_n_NICE_1000_n_TORAX_30[1])

rmse_n_NICE_1000_n_TORAX_50 = compute_rmse(phi_list_n_NICE_1000_n_TORAX_50[0], phi_list_n_NICE_1000_n_TORAX_50[1], rho_tor_norm_list_n_NICE_1000_n_TORAX_50[0], rho_tor_norm_list_n_NICE_1000_n_TORAX_50[1])

rmse_n_NICE_1000_n_TORAX_100 = compute_rmse(phi_list_n_NICE_1000_n_TORAX_100[0], phi_list_n_NICE_1000_n_TORAX_100[1], rho_tor_norm_list_n_NICE_1000_n_TORAX_100[0], rho_tor_norm_list_n_NICE_1000_n_TORAX_100[1])





rmse_n_NICE_50_list = [
    rmse_n_NICE_50_n_TORAX_10,
    rmse_n_NICE_50_n_TORAX_15,
    rmse_n_NICE_50_n_TORAX_20,
    rmse_n_NICE_50_n_TORAX_25,
    rmse_n_NICE_50_n_TORAX_30,
    rmse_n_NICE_50_n_TORAX_50,
    rmse_n_NICE_50_n_TORAX_100,
]


rmse_n_NICE_1000_list = [
    rmse_n_NICE_1000_n_TORAX_10,
    rmse_n_NICE_1000_n_TORAX_15,
    rmse_n_NICE_1000_n_TORAX_20,
    rmse_n_NICE_1000_n_TORAX_25,
    rmse_n_NICE_1000_n_TORAX_30,
    rmse_n_NICE_1000_n_TORAX_50,
    rmse_n_NICE_1000_n_TORAX_100,
]



n_TORAX_values = np.array([10, 15, 20, 25, 30, 50, 100])
n_NICE_values = np.array([25, 50, 100, 1000])


# Create 2D rmse_array, rows corresponding to n_NICE and column corresponding to n_TORAX

rmse_array = np.zeros((len(n_NICE_values), len(n_TORAX_values)))

rmse_array[1,:] = rmse_n_NICE_50_list
rmse_array[3,:] = rmse_n_NICE_1000_list
rmse_array[0,3] = rmse_n_NICE_25_n_TORAX_25
rmse_array[2,3] = rmse_n_NICE_100_n_TORAX_25

# Replace zeros with nan
rmse_array[rmse_array == 0] = np.nan


print(rmse_array)



# Plotting of heatmap

fig, ax = plt.subplots(figsize=(10, 6))
c = ax.imshow(rmse_array, cmap='YlGnBu', aspect='auto', interpolation='nearest')

# Add ticks
ax.set_xticks(np.arange(len(n_TORAX_values)))
ax.set_yticks(np.arange(len(n_NICE_values)))
ax.set_xticklabels(n_TORAX_values)
ax.set_yticklabels(n_NICE_values)
ax.set_xlabel("n_TORAX")
ax.set_ylabel("n_NICE")
ax.set_title("Consistency data through TORAX: RMSE between post_inverse and post_inverse $\phi$ for iteration 1, against NICE and TORAX profile resolution")

# Reverse y-axis
ax.invert_yaxis()

# Annotate values
for i in range(len(n_NICE_values)):
    for j in range(len(n_TORAX_values)):
        value = rmse_array[i, j]
        if not np.isnan(value):
            ax.text(j, i, f"{value:.2e}", ha="center", va="center", color="white" if value > np.nanmax(rmse_array)/2 else "black")

fig.colorbar(c, ax=ax, label="RMSE $\phi$ (Wb)")
plt.tight_layout()
plt.show()

