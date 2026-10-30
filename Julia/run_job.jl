# This function will be called by Python. It will take the following inputs:

# Filename: output filename
# Parameters: L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis
# Task: "thermalize" or "simulate"
# Number of sweeps
# One of:
#   (a) prev_file: filename to take states, operators from
#   (b) state, M: create an initial state with M padding identity operators M copies of state
# One of:
#   (a) use_prev_rng: true or false
#   (b) seed: create a new rng with this seed

function run_job(filename, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, task, num_sweeps; prev_file = "", state=0, M=0, use_prev_rng=false, seed=0)
    # Create / recover the configuration
    if prev_file == ""
        states = [state for _ in 1:M]
        operators = [[0, 0] for _ in 1:M]
    else
        h5open(prev_file, "r") do file
            states = read(file["configuration/states"])
            operators_mat = read(file["configuration/operators"])
            operators = [operators_mat[:, i] for i in axes(operators_mat, 2)]
        end
    end

    # Create / recover the RNG
    if use_prev_rng
        h5open(prev_file, "r") do file
            seed = read(file["metadata/seed"])
            rng_state = read(file["configuration/rng_state"])
            rng = deserialize(IOBuffer(rng_state))
        end
    else
        rng = Xoshiro(seed)
    end

    # Create / recover the diagonal operator list for this value of L
    # [REVISION 2026-10-23] Always build the list from the (corrected) lattice functions. A cache written by the
    # original code contains the old, wrong vertex masks and would silently bring the bug back, so a cached file is
    # kept only if it is identical to the freshly built list; otherwise it is overwritten with a warning.
    diagonals = find_all_diagonals(L)
    diagonals_filename = "Precomputed/diagonals_L$(L).h5"
    if isfile(diagonals_filename)
        cached = h5open(diagonals_filename, "r") do file
            read(file["diagonals"])
        end
        if cached != reduce(hcat, diagonals)
            @warn "Cached $(diagonals_filename) differs from the current lattice definition; overwriting it."
            rm(diagonals_filename)
        end
    end
    if !isfile(diagonals_filename)
        mkpath(dirname(diagonals_filename))
        h5open(diagonals_filename, "w") do file
            file["diagonals"] = reduce(hcat, diagonals)
        end
    end

    # Run the task
    if task == "thermalize"
        thermalize(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename)
    elseif task == "simulate"
        simulate(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename)
    elseif task == "test"
        test_offdiagonal(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename)
    else
        error("Not a valid task.")
    end
end