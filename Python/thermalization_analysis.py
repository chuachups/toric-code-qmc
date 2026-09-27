import h5py
import numpy as np
import matplotlib.pyplot as plt

# Recover the raw data
filename = "RawData/__L_2/__beta_1.000__J1_1.000__J2_1.000/__basis_X/run_0/thermalize_0.h5"
with h5py.File(filename, "r") as file:
    # Measurements
    M_arr = np.array(file["measurements/M"])
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

# Plot
sweeps = np.arange(0, num_sweeps)
plt.plot(sweeps, n_arr, label="n")
plt.plot(sweeps, M_arr, label="M") 

plt.xlabel("Sweep number") 
plt.ylabel("Operator count") 
plt.legend()
plt.show()
