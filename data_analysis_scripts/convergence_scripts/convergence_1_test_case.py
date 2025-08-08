import imas
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider
matplotlib.use("TkAgg")



""" The aim of this script is to be able to visualise convergence per iteration in terms of Root-Mean-Square-Error (RMSE) for a single test case """



def load_equilibria(number_of_iterations, test_case_number, base_path):
    """
    Load IMAS equilibrium entries for a specific test case across initial and multiple iteration steps.

    Args:
        number_of_iterations (int): Number of iterations to load.
        test_case_number (int): Test case number (j).
        base_path (str): Base path where test folders are stored.

    Returns:
        List of equilibrium objects, with the initial data at index 0.
    """
    equilibria = []

    # Construct base test case path
    test_case_path = os.path.join(base_path, f"test_{test_case_number}")

    # Load initial data
    initial_path = os.path.join(test_case_path, "initial_data")
    initial_uri = f"imas:hdf5?path={initial_path}"
    initial_entry = imas.DBEntry(initial_uri, "r")
    initial_eq = initial_entry.get("equilibrium", lazy=True)
    equilibria.append(initial_eq)

    # Load iteration data
    for i in range(1, number_of_iterations + 1):
        iteration_path = os.path.join(test_case_path, f"iteration_{i}", "post_torax")
        uri = f"imas:hdf5?path={iteration_path}"
        entry = imas.DBEntry(uri, "r")
        eq = entry.get("equilibrium", lazy=True)
        equilibria.append(eq)

    return equilibria



def compute_rmse(array1, array2, x1=None, x2=None):
    """
    Computes RMSE between two 1D profiles, interpolating to their overlapping domain if needed.

    If the profiles are already matching in shape, RMSE is directly computed.

    If the profiles differ in shape, the denser is interpolated onto the coarser within their overlapping range.

    Args:
        array1 (np.ndarray): First profile (e.g. psi_1).
        array2 (np.ndarray): Second profile (e.g. psi_2).
        x1 (np.ndarray, optional): Grid for array1 (e.g. rho_1).
        x2 (np.ndarray, optional): Grid for array2 (e.g. rho_2).

    Returns:
        float: RMSE within overlapping domain.

    Raises:
        ValueError: If no overlapping domain is found or if grid information is missing when needed.
    """
    array1 = np.asarray(array1)
    array2 = np.asarray(array2)

    # If profiles are directly matching in shape
    if array1.shape == array2.shape:
        return np.sqrt(np.mean((array1 - array2) ** 2))

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
    interpolator1 = interp1d(x1, array1, kind='linear')
    interpolator2 = interp1d(x2, array2, kind='linear')

    array1_new = interpolator1(new_grid)
    array2_new = interpolator2(new_grid)

    # Calculate RMSE on overlapping domain
    rmse = np.sqrt(np.mean((array1_new - array2_new) ** 2))

    return rmse



equilibria_test_15 = load_equilibria(number_of_iterations = 10, test_case_number = 13, base_path = "/home/ITER/vanschr/investigation_nice_inv_torax_coupling_scheme")

psi_list_15 = [
    eq.time_slice[0].profiles_1d.psi for eq in equilibria_test_15
]

rho_tor_norm_list_15 = [
    eq.time_slice[0].profiles_1d.rho_tor_norm for eq in equilibria_test_15
]


rmse_psi_list_15 = [
    compute_rmse(psi_list_15[i-1], psi_list_15[i],
                 rho_tor_norm_list_15[i-1], rho_tor_norm_list_15[i])
    for i in range(1, len(psi_list_15))
]


# X-axis: iteration numbers (1-based index for comparison i vs i-1)
iteration_numbers = list(range(1, len(rmse_psi_list_15)+1))

plt.figure(figsize=(8, 5))
plt.scatter(iteration_numbers[1:], rmse_psi_list_15[1:], color='blue', marker='o')
plt.plot(iteration_numbers[1:], rmse_psi_list_15[1:], linestyle='--', color='gray', alpha=0.7)

plt.title("RMSE of $\psi$ profiles between successive iterations")
plt.xlabel("Iteration")
plt.ylabel("RMSE of $\psi$")
#plt.grid(True)
plt.tight_layout()
plt.show()