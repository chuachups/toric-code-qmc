# This Python function will take the simple job description (stored in jobs.py) and construct the complete command required by Julia
#
# Each job specifies the following information
# 1. Parameters: L, beta, J_1, J_2, basis
# 2. Run number
# 3. Task: "thermalize" or "simulate"
# 4. Number of sweeps
# 5. Optional for thermalization warm start: L, beta, J_1, J_2, basis, run number of the previous run to use as a warm start
#
# The filename convention is, for example, "RawData/__L_2/__beta_1__J1_1__J2_1/__basis_Z/run_0/thermalize_0"
# The logic is:
# # New run:
#   - choose a new RNG seed
#   - if a warm-start source is given, use its latest simulation, or latest thermalization if no simulation exists
#   - otherwise create a random initial state and determine M based on L, beta
# Existing run:
#   - thermalize continues the latest thermalization
#   - simulate continues the latest simulation, or the latest thermalization if no simulation exists

import os
import random
import math

def make_job(job, pbs=False):
    # Unpack
    L = job["L"]
    beta = job["beta"]
    J_1 = job["J_1"]
    J_2 = job["J_2"]
    basis = job["basis"]
    run = job["run"]
    task = job["task"]
    num_sweeps = job["num_sweeps"]

    # Automatically determine the constant shifts
    C_1 = 2
    C_2 = 2
    C_3 = 1
    if (J_1 == 0 and basis == "Z") or (J_2 == 0 and basis == "X"):
        C_4 = 1
    else:
        C_4 = 0

    # Construct filename and determine the latest entry in the given run
    filename = make_file_name(L, beta, J_1, J_2, basis, run, task)
    entry_num = latest_entry(filename)

    # Construct the full job description, which is a list of the following parameters

    # For "thermalize"
    if task == "thermalize":
        # Check that no simulate exists
        if latest_entry(make_file_name(L, beta, J_1, J_2, basis, run, "simulate")) != -1:
            raise ValueError(f"Cannot thermalize run, simulation data already exists.")

        # If no previous thermalization exists, create a random rng seed
        if entry_num == -1:
            rng_seed = random.randrange(2**63)

            # If a warm start path is given, use it
            if "previous" in job:
                # Unpack
                prev = job["previous"]
                L_prev = prev["L"]
                beta_prev = prev["beta"]
                J_1_prev = prev["J_1"]
                J_2_prev = prev["J_2"]
                basis_prev = prev["basis"]
                run_prev = prev["run"]

                # Find the latest simulation
                filename_prev = make_file_name(L_prev, beta_prev, J_1_prev, J_2_prev, basis_prev, run_prev, "simulate")
                entry_num_prev = latest_entry(filename_prev)

                # If it doesn't exist, find the latest thermalization
                if entry_num_prev == -1:
                    # Find the latest thermalization
                    filename_prev = make_file_name(L_prev, beta_prev, J_1_prev, J_2_prev, basis_prev, run_prev, "thermalize")
                    entry_num_prev = latest_entry(filename_prev)

                    # If it doesn't exist, raise an error
                    if entry_num_prev == -1:
                        raise FileNotFoundError("No existing entries for the specified warm start run yet.")

                output_file = filename + str(entry_num + 1) + ".h5"
                prev_file = filename_prev + str(entry_num_prev) + ".h5"
                state = 0
                M = 0
                use_prev_rng = False
                seed = rng_seed

            # If no warm start is given, create a random initial state and determine M
            else:
                state_init = random.randrange(2**(2 * L**2))
                M_init = math.ceil(beta * L**2)

                output_file = filename + str(entry_num + 1) + ".h5"
                prev_file = ""
                state = state_init
                M = M_init
                use_prev_rng = False
                seed = rng_seed

        # If a previous thermalization exists, use it
        else:
            # If a warm start path is specified, raise an error
            if "previous" in job:
                raise ValueError("Cannot use warm start path because data already exists.")

            output_file = filename + str(entry_num + 1) + ".h5"
            prev_file = filename + str(entry_num) + ".h5"
            state = 0
            M = 0
            use_prev_rng = True
            seed = 0
    
    # If the task is "simulate", determine the latest simulation
    elif task == "simulate":
        # If a warm start path is specified, raise an error
        if "previous" in job:
            raise ValueError("Cannot use warm start path for simulate.")
            
        # If no previous simulation exists, use the latest thermalization
        if entry_num == -1:
            # Find the latest thermalization
            filename_prev = make_file_name(L, beta, J_1, J_2, basis, run, "thermalize")
            entry_num_prev = latest_entry(filename_prev)

            # If it doesn't exist, raise an error
            if entry_num_prev == -1:
                raise FileNotFoundError("No existing thermalization for this set of parameters yet.")

            output_file = filename + str(entry_num + 1) + ".h5"
            prev_file = filename_prev + str(entry_num_prev) + ".h5"
            state = 0
            M = 0
            use_prev_rng = True
            seed = 0

        # If a previous simulation exists, use it
        else:
            output_file = filename + str(entry_num + 1) + ".h5"
            prev_file = filename + str(entry_num) + ".h5"
            state = 0
            M = 0
            use_prev_rng = True
            seed = 0

    else:
        raise ValueError(f"Invalid task: {task}")

    # Construct the full command
    if pbs == False:
        command = [
                "julia",
                "Julia/main.jl",
        
                "--filename", output_file,
                "--L", str(L),
                "--beta", str(beta),
                "--J1", str(J_1),
                "--J2", str(J_2),
                "--C1", str(C_1),
                "--C2", str(C_2),
                "--C3", str(C_3),
                "--C4", str(C_4),
                "--basis", basis,
                "--task", task,
                "--num-sweeps", str(num_sweeps),
            ]
    else:
        command = [
                "julia",
                "--project=.",
                "Julia/main.jl",
        
                "--filename", output_file,
                "--L", str(L),
                "--beta", str(beta),
                "--J1", str(J_1),
                "--J2", str(J_2),
                "--C1", str(C_1),
                "--C2", str(C_2),
                "--C3", str(C_3),
                "--C4", str(C_4),
                "--basis", basis,
                "--task", task,
                "--num-sweeps", str(num_sweeps),
            ]
    if prev_file != "":
        command += ["--prev-file", prev_file]
    else:
        command += ["--state", str(state), "--M", str(M)]
    if use_prev_rng:
        command += ["--use-prev-rng"]
    else:
        command += ["--seed", str(seed)]

    return command

# Helper function, create the filename
def make_file_name(L, beta, J_1, J_2, basis, run, task):
    return f"RawData/__L_{L}/__beta_{beta:.3f}__J1_{J_1:.3f}__J2_{J_2:.3f}/__basis_{basis}/run_{run}/{task}_"

# Helper function, determine the latest entry of a run
def latest_entry(filename_prefix):
    entry_num = 0
    while os.path.isfile(filename_prefix + str(entry_num) + ".h5"):
        entry_num += 1
    return entry_num - 1
    