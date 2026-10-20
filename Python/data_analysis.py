import h5py
import numpy as np

# This function computes the desired observable given part of the data
def calculate_observable(observable, data, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis):
    n_arr, n2_arr, n_A_arr, n_B_arr, n_zz_arr, n_xx_arr, m_arr, m_2_arr, m_4_arr, chi_M_arr, n_j_arr, n2_j_arr, S_kmin_arr = data
    if observable == "E":
        if basis == "Z":
            return -np.mean(n_arr) / beta + (C_1 + C_3 + 4 * J_1 * C_2 + 4 * C_4) * L**2
        else:
            return -np.mean(n_arr) / beta + (C_1 + C_3 + 4 * J_2 * C_2 + 4 * C_4) * L**2
        # return -(np.mean(n_A_arr) + np.mean(n_B_arr) + np.mean(n_zz_arr) + np.mean(n_xx_arr)) / beta + (C_1 + 4 * J_1 * C_2) * L**2
    elif observable == "M^2":
        return np.mean(m_2_arr)
    elif observable == "M^4":
        return np.mean(m_4_arr)
    elif observable == "C":
        return np.mean(n2_arr) - np.mean(n_arr)**2 - np.mean(n_arr)
    elif observable == "chi_M":
        return np.mean(chi_M_arr)
    elif observable == "chi_J":
        return (np.mean(n2_j_arr) - np.mean(n_j_arr)**2 - np.mean(n_j_arr)) / (beta * 2 * L**2 * J_1**2)
    elif observable == "xi_M":
        S_0 = np.mean(m_2_arr) * 2*L**2
        if S_0 > np.mean(S_kmin_arr):
            return (1 / (2 * np.sin(np.pi / L))) * np.sqrt(S_0 / np.mean(S_kmin_arr) - 1)
        else:
            return 0
    elif observable == "Binder":
        return 1 - np.mean(m_4_arr) / (3 * np.mean(m_2_arr)**2)
    elif observable == "A":
        if basis == "Z":
            return np.mean(n_A_arr) / (beta * L**2)
        else:
            return np.mean(n_A_arr) / (beta * L**2) - C_1
    elif observable == "B":
        if basis == "Z":
            return np.mean(n_B_arr) / (beta * L**2) - C_1
        else:
            return np.mean(n_B_arr) / (beta * L**2)
    elif observable == "ZZ":
        if basis == "Z":
            return np.mean(n_zz_arr) / (beta * J_1 * 4 * L**2) - C_2
        else:
            return np.mean(n_zz_arr) / (beta * J_1 * 4 * L**2)
    elif observable == "XX":
        if J_2 == 0:
            return 0
        elif basis == "Z":
            return np.mean(n_xx_arr) / (beta * J_2 * 4 * L**2)
        else:
            return np.mean(n_xx_arr) / (beta * J_2 * 4 * L**2) - C_2

# This function performs Jackknife given the desired observable and all the data
def jackknife(observable, data, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis):
    # Initially take b to be 5
    # Then, double it until the estimated error converges
    b = 5

    # Do a machimum of 30 doublings
    err = 0
    for i in range(30):
        # Separate the data into B blocks of size b, and compute the means with each block left out
        data_wo_blk_means, num_blks = wo_blk_means(data, b)

        # Stop if too few blocks
        if num_blks < 10:
            break

        # Calculate the observable
        results = np.zeros(num_blks)
        for i in range(num_blks):
            results[i] = calculate_observable(observable, [qty[i] for qty in data_wo_blk_means], L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis)

        # Compute mean and error
        mean = np.mean(results)
        err_new = np.sqrt(num_blks - 1) * np.std(results)

        # Check for convergence
        if err_new == 0:
            break
        if abs((err_new - err) / err_new) < 5 * 10**(-2):
            err = err_new
            break
        
        # Double b and update U_err
        err = err_new
        b *= 2

    return calculate_observable(observable, data, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis), err

# This function separates the data into B blocks of size b and computes the means with each block left out
def wo_blk_means(data, b):
    data_wo_blk_means = []
    for qty in data:
        # Remove first few entries to ensure that B is an integer.
        qty = qty[qty.size % b:]
        # Then resize the array. Each row corresponds to a block.
        qty_blk = qty.reshape(-1, b)
        # Then take the means across each row.
        qty_blk_means = np.mean(qty_blk, axis=1)
        # Compute the means with each block left out
        num_blks = qty_blk_means.size
        qty_wo_blk_means = (num_blks * np.mean(qty_blk_means) - qty_blk_means) / (num_blks - 1)
        # Append
        data_wo_blk_means.append(qty_wo_blk_means)

    return data_wo_blk_means, num_blks

# This function formats the result, rounding to the error
def format_result(mean, err):
    if err == 0:
        return f"{mean} ± {0}"
    else:
        decimal_place = int(np.floor(np.log10(abs(err)))) - 1

        err_round = round(err, -decimal_place)
        mean_round = round(mean, -decimal_place)

        return f"{mean_round} ± {err_round}"

# This function computes all the observables with error bars given the filename, and returns a dictionary
def compute_all(filename):
    with h5py.File(filename, "r") as file:
        # Measurements
        M = file["measurements/M"][()]
        n_arr = np.array(file["measurements/n"])
        n_A_arr = np.array(file["measurements/n_A"])
        n_B_arr = np.array(file["measurements/n_B"])
        n_zz_arr = np.array(file["measurements/n_zz"])
        n_xx_arr = np.array(file["measurements/n_xx"])
        m_arr = np.array(file["measurements/m"])
        m_2_arr = np.array(file["measurements/m_2"])
        m_4_arr = np.array(file["measurements/m_4"])
        S_kmin_arr = np.array(file["measurements/S_kmin"])

        # Parameters
        L = file["parameters/L"][()]
        beta = file["parameters/beta"][()]
        J_1 = file["parameters/J_1"][()]
        J_2 = file["parameters/J_2"][()]
        C_1 = file["parameters/C_1"][()]
        C_2 = file["parameters/C_2"][()]
        C_3 = file["parameters/C_3"][()]
        C_4 = file["parameters/C_4"][()]

        # Metadata
        basis = file["metadata/basis"][()].decode()
        num_sweeps = file["metadata/num_sweeps"][()]

    # Pre-compute
    n2_arr = n_arr**2
    chi_M_arr = 2 * L**2 * beta * (n_arr * m_arr**2 + m_2_arr) / ((n_arr + 1))
    n_j_arr = n_zz_arr + n_xx_arr
    n2_j_arr = n_j_arr**2

    data = [n_arr, n2_arr, n_A_arr, n_B_arr, n_zz_arr, n_xx_arr, m_arr, m_2_arr, m_4_arr, chi_M_arr, n_j_arr, n2_j_arr, S_kmin_arr]

    means = {}
    errs = {}
    for observable in ["E", "M^2", "M^4", "C", "chi_M", "chi_J", "xi_M", "Binder", "A", "B", "ZZ", "XX"]:
        mean, err = jackknife(observable, data, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis)
        means[observable] = mean
        errs[observable] = err

    return means, errs

# This function computes all the observables given the filename, and prints the results
def print_all(filename):
    with h5py.File(filename, "r") as file:
        # Measurements
        M = file["measurements/M"][()]
        n_arr = np.array(file["measurements/n"])
        n_A_arr = np.array(file["measurements/n_A"])
        n_B_arr = np.array(file["measurements/n_B"])
        n_zz_arr = np.array(file["measurements/n_zz"])
        n_xx_arr = np.array(file["measurements/n_xx"])
        m_arr = np.array(file["measurements/m"])
        m_2_arr = np.array(file["measurements/m_2"])
        m_4_arr = np.array(file["measurements/m_4"])
        S_kmin_arr = np.array(file["measurements/S_kmin"])

        # Parameters
        L = file["parameters/L"][()]
        beta = file["parameters/beta"][()]
        J_1 = file["parameters/J_1"][()]
        J_2 = file["parameters/J_2"][()]
        C_1 = file["parameters/C_1"][()]
        C_2 = file["parameters/C_2"][()]
        C_3 = file["parameters/C_3"][()]
        C_4 = file["parameters/C_4"][()]

        # Metadata
        basis = file["metadata/basis"][()].decode()
        num_sweeps = file["metadata/num_sweeps"][()]

    # Pre-compute
    n2_arr = n_arr**2
    chi_M_arr = 2 * L**2 * beta * (n_arr * m_arr**2 + m_2_arr) / ((n_arr + 1))
    n_j_arr = n_zz_arr + n_xx_arr
    n2_j_arr = n_j_arr**2

    data = [n_arr, n2_arr, n_A_arr, n_B_arr, n_zz_arr, n_xx_arr, m_arr, m_2_arr, m_4_arr, chi_M_arr, n_j_arr, n2_j_arr, S_kmin_arr]

    # Compute observables
    print(f"L={L}, beta={beta}, J_1={J_1}, J_2={J_2}, basis={basis}")
    print("Cutoff:", M)
    for observable in ["E", "M^2", "M^4", "C", "chi_M", "chi_J", "xi_M", "Binder", "A", "B", "ZZ", "XX"]:
        mean, err = jackknife(observable, data, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis)
        print(observable, ":", format_result(mean, err))

print_all("RawData/__L_2/__beta_5.000__J1_0.100__J2_0.100/__basis_Z/run_1/simulate_1.h5")