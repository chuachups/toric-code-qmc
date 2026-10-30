# [REVISION 2026-10-23] New keyword arguments m_n_arr, m_2_n_arr: <m_z> and <m_z^2> over the n propagated states,
# used by data_analysis.py for the chi_M estimator. They are written only if non-empty.
function write_output_file(states_init, operators_init, states, operators, M_data, n_arr, n_A_arr, n_B_arr, n_zz_arr, n_xx_arr, m_arr, m_2_arr, m_4_arr, S_kmin_arr, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename; m_n_arr=Float64[], m_2_n_arr=Float64[])
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
        file["measurements/S_kmin"] = S_kmin_arr
        # [REVISION 2026-10-23] Averages over the n propagated states (chi_M estimator)
        if !isempty(m_n_arr)
            file["measurements/m_n"] = m_n_arr
            file["measurements/m_2_n"] = m_2_n_arr
        end

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
        #file["metadata/version"] = strip(read(`git rev-parse HEAD`, String))
        file["metadata/date"] = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSxxx")
        file["metadata/basis"] = basis
        file["metadata/num_sweeps"] = num_sweeps
        file["metadata/seed"] = seed
        file["metadata/states_init"] = states_init
        file["metadata/operators_init"] = reduce(hcat, operators_init)
    end
end