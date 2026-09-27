import h5py
import numpy as np

# Recover the raw data
filename = "RawData/__L_2/__beta_1.000__J1_1.000__J2_1.000/__basis_X/run_0/simulate_0.h5"
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
    S_0_arr = np.array(file["measurements/S_0"])
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
    basis = file["metadata/basis"][()]
    num_sweeps = file["metadata/num_sweeps"][()]

# Compute observables
print("Cutoff:", M)

print("Energy:", -np.mean(n_arr) / beta + (C_1 + C_3 + 4 * J_1 * C_2 + 4 * C_4) * L**2)
print("Energy:", -(np.mean(n_A_arr) + np.mean(n_B_arr) + np.mean(n_zz_arr) + np.mean(n_xx_arr)) / beta + (C_1 + 4 * J_1 * C_2) * L**2)

print("M_z^2:", np.mean(m_2_arr))
print("M_z^4:", np.mean(m_4_arr))

print("Heat capacity:", np.mean(n_arr**2) - np.mean(n_arr)**2 - np.mean(n_arr))

xi_z_arr = 2 * L**2 * beta * ((n_arr * m_arr)**2 + n_arr * m_2_arr) / (n_arr * (n_arr + 1))
print("Susceptibility (Z):", np.mean(xi_z_arr))

n_j_arr = n_zz_arr + n_xx_arr
print("Susceptibility (J):", (np.mean(n_j_arr**2) - np.mean(n_j_arr)**2 - np.mean(n_j_arr)) / (beta * 2 * L**2 * J_1**2))

print("Correlation length:", (1 / (2 * np.sin(np.pi / L))) * np.sqrt(np.mean(S_0_arr) / np.mean(S_kmin_arr) - 1))

if basis == "Z":
    print("A:", np.mean(n_A_arr) / (beta * L**2))
    print("B:", np.mean(n_B_arr) / (beta * L**2) - C_1)
    print("ZZ:", np.mean(n_zz_arr) / (beta * J_1 * 4 * L**2) - C_2)
    print("XX:", np.mean(n_xx_arr) / (beta * J_2 * 4 * L**2))
else:
    print("A:", np.mean(n_A_arr) / (beta * L**2) - C_1)
    print("B:", np.mean(n_B_arr) / (beta * L**2))
    print("ZZ:", np.mean(n_zz_arr) / (beta * J_1 * 4 * L**2))
    print("XX:", np.mean(n_xx_arr) / (beta * J_2 * 4 * L**2) - C_2)