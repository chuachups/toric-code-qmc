import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Generate a plot for both n and M^2
def gen_all_plots():
    for var in ["n", "m2"]:
        # Loop through all files in RawData
        raw_dir = Path("RawData")
        analysis_dir = Path("Analysis")
        files = sorted(raw_dir.rglob("thermalize_*.h5"))
        for filename in files:
            relative_path = filename.relative_to(raw_dir)
            output_path = analysis_dir / relative_path.with_stem(relative_path.stem + f"_{var}").with_suffix(".png")

            # Skip if the plot already exists
            if output_path.exists():
                continue

            # Make sure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Generate the plot
            with h5py.File(filename, "r") as file:
                # Get raw data
                if var == "n":
                    M_arr = np.array(file["measurements/M"])
                    arr = np.array(file["measurements/n"])
                else:
                    arr = np.array(file["measurements/m_2"])
                L = file["parameters/L"][()]
                beta = file["parameters/beta"][()]
                J_1 = file["parameters/J_1"][()]
                J_2 = file["parameters/J_2"][()]
                basis = file["metadata/basis"][()].decode()
                num_sweeps = file["metadata/num_sweeps"][()]

                # Compute the cumulative average
                sweeps = np.arange(1, num_sweeps + 1)
                cum_avg = np.cumsum(arr) / sweeps

                # Compute the window average
                window = 100
                win_avg = np.convolve(arr, np.ones(window) / window, mode="valid")
                win_sweeps = np.arange(window, num_sweeps + 1)

                # Plot
                plt.figure()
                if var == "n":
                    plt.plot(sweeps, arr, label="Instantaneous n")
                    plt.plot(win_sweeps, win_avg, label="Window n")
                    plt.plot(sweeps, cum_avg, label="Cumulative n")
                    plt.plot(sweeps, M_arr, label="M") 
                    plt.ylabel("Operator count") 
                else:
                    plt.plot(sweeps, arr, label="Instantaneous M^2")
                    plt.plot(win_sweeps, win_avg, label="Window M^2")
                    plt.plot(sweeps, cum_avg, label="Cumulative M^2")
                    plt.ylabel("Magnetization squared") 
                plt.xlabel("Sweep number") 
                plt.title(f"L={L}, beta={beta}, J_1={J_1}, J_2={J_2}, basis={basis}")
                plt.legend()    

                plt.savefig(output_path, dpi=150, bbox_inches="tight")
                plt.close()

                print(f"Generated {output_path}")

# Custom: Generate a plot for a specific file
def gen_plot(filename, var):
    # Generate the plot
    with h5py.File(filename, "r") as file:
        # Get raw data
        if var == "n":
            M_arr = np.array(file["measurements/M"])
            arr = np.array(file["measurements/n"])
        else:
            arr = np.array(file["measurements/m_2"])
        L = file["parameters/L"][()]
        beta = file["parameters/beta"][()]
        J_1 = file["parameters/J_1"][()]
        J_2 = file["parameters/J_2"][()]
        basis = file["metadata/basis"][()].decode()
        num_sweeps = file["metadata/num_sweeps"][()]

        # Compute the cumulative average
        sweeps = np.arange(1, num_sweeps + 1)
        cum_avg = np.cumsum(arr) / sweeps

        # Compute the window average
        window = 100
        win_avg = np.convolve(arr, np.ones(window) / window, mode="valid")
        win_sweeps = np.arange(window, num_sweeps + 1)

        # Plot
        if var == "n":
            plt.plot(sweeps, arr, label="Instantaneous n")
            plt.plot(win_sweeps, win_avg, label="Window n")
            plt.plot(sweeps, cum_avg, label="Cumulative n")
            plt.plot(sweeps, M_arr, label="M") 
            plt.ylabel("Operator count") 
        else:
            plt.plot(sweeps, arr, label="Instantaneous M^2")
            plt.plot(win_sweeps, win_avg, label="Window M^2")
            plt.plot(sweeps, cum_avg, label="Cumulative M^2")
            plt.ylabel("Magnetization squared") 
        plt.xlabel("Sweep number") 
        plt.title(f"L={L}, beta={beta}, J_1={J_1}, J_2={J_2}, basis={basis}")
        plt.legend()

        plt.show()

#gen_all_plots()
gen_plot("RawData/__L_3/__beta_5.000__J1_0.100__J2_0.000/__basis_Z/run_0/simulate_1.h5", "m2")