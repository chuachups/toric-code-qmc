function write_output_file(states_init, operators_init, states, operators, M_data, n_arr, n_A_arr, n_B_arr, n_zz_arr, n_xx_arr, m_arr, m_2_arr, m_4_arr, S_0_arr, S_kmin_arr, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename)
    # Convert the RNG state to bytes
    io = IOBuffer()
    serialize(io, rng)
    rng_state = take!(io)

    mkpath(dirname(filename))
    h5open(filename, "w") do file
        # Store the current configuration (needed to resume the run for warm starts)
        file["configuration/states"] = states
        file["configuration/operators"] = reduce(hcat, operators)
        file["configuration/rng_state"] = rng_state

        # Store the raw data
        file["measurements/M"] = M_data
        file["measurements/n"] = n_arr
        file["measurements/n_A"] = n_A_arr
        file["measurements/n_B"] = n_B_arr
        file["measurements/n_zz"] = n_zz_arr
        file["measurements/n_xx"] = n_xx_arr
        file["measurements/m"] = m_arr
        file["measurements/m_2"] = m_2_arr
        file["measurements/m_4"] = m_4_arr
        file["measurements/S_0"] = S_0_arr
        file["measurements/S_kmin"] = S_kmin_arr

        # Store the parameters
        file["parameters/L"] = L
        file["parameters/beta"] = beta
        file["parameters/J_1"] = J_1
        file["parameters/J_2"] = J_2
        file["parameters/C_1"] = C_1
        file["parameters/C_2"] = C_2
        file["parameters/C_3"] = C_3
        file["parameters/C_4"] = C_4

        # Store metadata for reproducibility
        file["metadata/version"] = strip(read(`git rev-parse HEAD`, String))
        file["metadata/date"] = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSxxx")
        file["metadata/basis"] = basis
        file["metadata/num_sweeps"] = num_sweeps
        file["metadata/seed"] = seed
        file["metadata/states_init"] = states_init
        file["metadata/operators_init"] = reduce(hcat, operators_init)
    end
end