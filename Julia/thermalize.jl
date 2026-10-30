# Given the configuration, perform a given number of thermalization sweeps
function thermalize(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename)
    # Keep a copy of the initial configuration
    states_init = copy(states)
    operators_init = deepcopy(operators)
    
    # If the basis is X, swap J_1, J_2 and (later) return in the correct order
    if basis == "X"
        J_1, J_2 = J_2, J_1
    end

    # This is essentially the same as simulate(), except we allow a variable M
    M_arr = Float64[]
    n_arr = Float64[]
    n_A_arr = Float64[]
    n_B_arr = Float64[]
    n_zz_arr = Float64[]
    n_xx_arr = Float64[]
    mz_arr = Float64[]
    mz_2_arr = Float64[]
    mz_4_arr = Float64[]
    S_kmin_arr = Float64[]
    m_n_arr = Float64[]     # [REVISION 2026-10-23] <m_z> over the n propagated states (for chi_M)
    m_2_n_arr = Float64[]   # [REVISION 2026-10-23] <m_z^2> over the n propagated states (for chi_M)

    # If J_2 = 0, the initial configuration should have no XX operators
    if J_2 == 0 && count(op -> op[1] == 4, operators) != 0
        error("Error: Initial configuration must not contain any XX-type operators for J_2 = 0.")
    end

    for i in 1:num_sweeps
        # Perform sweep
        states, operators = qmc_sweep(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, rng, false)
        # Perform measurements
        n, n_A, n_B, n_zz, n_xx, mz, mz_2, mz_4, S_kmin = measure_all(states, operators, L)
        push!(M_arr, length(states))
        push!(n_arr, n)
        push!(n_A_arr, n_A)
        push!(n_B_arr, n_B)
        if J_1 == 0
            push!(n_zz_arr, 0)
        else
            push!(n_zz_arr, n_zz)
        end
        push!(n_xx_arr, n_xx)
        push!(mz_arr, mz)
        push!(mz_2_arr, mz_2)
        push!(mz_4_arr, mz_4)
        push!(S_kmin_arr, S_kmin)
        # [REVISION 2026-10-23] Averages over the n propagated states, for the chi_M estimator
        m_n, m_2_n = measure_m_propagated(states, operators, L)
        push!(m_n_arr, m_n)
        push!(m_2_n_arr, m_2_n)
        # Update the cutoff
        # [REVISION 2026-10-23] Two fixes: (i) padding operators are [0, 0] (was [0], which made
        # reduce(hcat, operators) in write_output_file throw a DimensionMismatch); (ii) the new slots come after the
        # last operator, so their state is the state AFTER operator M, which by periodicity is states[1] (was
        # states[end], the state BEFORE operator M, which is wrong whenever operator M is XXXX or XX).
        state_after_last = states[1]
        for _ in 1:max(0, ceil(Int, 4/3 * n - length(states)))
            push!(states, state_after_last)
            push!(operators, [0, 0])
        end
    end

    if basis == "X"
        write_output_file(states_init, operators_init, states, operators, M_arr, n_arr, n_B_arr, n_A_arr, n_xx_arr, n_zz_arr, mz_arr, mz_2_arr, mz_4_arr, S_kmin_arr, L, beta, J_2, J_1, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename; m_n_arr=m_n_arr, m_2_n_arr=m_2_n_arr)
    else
        write_output_file(states_init, operators_init, states, operators, M_arr, n_arr, n_A_arr, n_B_arr, n_zz_arr, n_xx_arr, mz_arr, mz_2_arr, mz_4_arr, S_kmin_arr, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, basis, num_sweeps, rng, seed, filename; m_n_arr=m_n_arr, m_2_n_arr=m_2_n_arr)
    end
end