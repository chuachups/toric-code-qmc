# Brief Summary
###############################################################################################
# This file contains all the functions required to run the Quantum Monte Carlo (QMC) algorithm.
# The QMC used here is based on a stochastic series expansion (SSE) of the partition function in
# the computational basis. Monte Carlo is used to sample the terms of the expansion, which will
# be referred to as "configurations". Each configuration is specified by n computational basis
# states |alpha(p)> and n operators H_p  such that |alpha(p+1)> = H_p |alpha(p)>. Each operator 
# can be any term in the Hamiltonian, and probability weight of a particular configuration is 
# given by (beta^n / n!) * <alpha(0)|H_n...H_1|alpha(0)>. In practice, there will be a fixed 
# cutoff M states and operators, and variable n is achieved by allowing H_p to be an identity
# operator to serve as padding.

# The original Hamiltonian is:
#
# -H = sum_v (XXXX) + sum_p (ZZZZ) + J_1 sum_b (ZZ) + J_2 sum_b (XX)
#
# To facilitate the algorithm, several constant terms are added. The full Hamiltonian is:
#
# -H = sum_v (XXXX) + sum_v (C_3 IIII) + sum_p (ZZZZ + C_1) + sum_b (J_1 (ZZ + C_2)) + sum_b (J_2 XX)
#
# The operators in the configuration can be any of the vertex terms XXXX or IIII, plaquette terms 
# ZZZZ, bond terms ZZ or XX, or the padding identity I. 
# The diagonal updates attempt to replace each padding identity I with a random diagonal operator 
# ZZZZ, IIII, or ZZ, and vice versa. (I <-> ZZZZ, IIII, ZZ)
# The off-diagonal updates attempt to replace diagonal operators with the corresponding off-diagonal
# operators (IIII <-> XXXX, ZZ <-> XX). This is done simultaneously for all operators linked together
# by the spins they act on (referred to as a cluster). The spins in the cluster are flipped (up <-> down).

# Conventions:
# 1. Each spin is assigned an index i from 0 to 2L^2 - 1. 
# - For spins on horizontal links: i(x,y) = (x-1) + L(y-1) where (x,y) is the vertex to the left of 
#   the link.
# - For spins on vertical links: i(x,y) = L^2 + (x-1) + L(y-1) where (x,y) is the vertex above the
#   link.
# 2. Each computational basis state is assigned an integer from 0 to 2^(2L^2) - 1. 
# - The binary representation of this integer determines whether each spin is in the up or down
#   state. If bit i is 0 (1), then spin i is in the spin up (down) state.
# 3. Each operator is assigned a list of two integers. 
# - The first integer determines the type of the operator. The convention is:
#   0: I (padding)
#   1: ZZZZ (plaquette)
#   2: XXXX (vertex)
#   3: ZZ (bond)    (If J_1 = 0, then this is II)
#   4: XX (bond)    (If J_2 = 0, then this is absent)
#   5: IIII (vertex)
# - The binary representation of the second integer determines which spins the operator acts on.
#   If bit i is 1 (0), then the operator acts (does not act) on spin i.
###############################################################################################

# The code is organized into four main parts:
# 1. Several helper functions to generate the terms in the Hamiltonian
# 2. One main function to run the QMC SSE sampling + two helper functions
# 3. One main function to measure the raw data over many sweeps + several helper functions
# 4. Thermalization code which automatically determines M
# Data analysis is in Python.

# There are three kinds of terms in the Hamiltonian: vertices, plaquettes, and bonds. The (second)
# integers assigned to each are generated using the following three functions:

# Given the size L, return a list of all vertex operators.
function find_all_vertices(L)
    vertices = Int[]
    # Loop through all the vertices (x,y)
    for x in 1:L
        for y in 1:L
            vertex = 0
            # Find the top spin
            vertex |= 1 << ((x - 1) + L * mod(y - 2, L))
            # Find the bottom spin
            vertex |= 1 << ((x - 1) + L * (y - 1))
            # Find the left spin
            vertex |= 1 << (L^2 + mod(x - 2, L) + L * (y - 1))
            # Find the right spin
            vertex |= 1 << (L^2 + (x - 1) + L * (y - 1))
            push!(vertices, vertex)
        end
    end
    return vertices
end

# Given the size L, return a list of all plaquette operators.
function find_all_plaquettes(L)
    plaquettes = Int[]
    # Loop through all the plaquettes (x,y)
    for x in 1:L
        for y in 1:L
            plaquette = 0
            # Find the top spin
            plaquette |= 1 << ((x - 1) + L * (y - 1))
            # Find the bottom spin
            plaquette |= 1 << ((x - 1) + L * mod(y, L))
            # Find the left spin
            plaquette |= 1 << (L^2 + (x - 1) + L * (y - 1))
            # Find the right spin
            plaquette |= 1 << (L^2 + mod(x, L) + L * (y - 1))
            push!(plaquettes, plaquette)
        end
    end
    return plaquettes
end

# Given the size L, return a list of all the bond operators.
function find_all_bonds(L)
    bonds = Int[]
    # Loop through all the horizontal links
    for x in 1:L
        for y in 1:L
            # Add bond with top-left spin
            bond = 0
            bond |= 1 << ((x - 1) + L * (y - 1))
            bond |= 1 << (L^2 + (x - 1) + L * mod(y - 2, L))
            push!(bonds, bond)
            # Add bond with the bottom-left spin
            bond = 0
            bond |= 1 << ((x - 1) + L * (y - 1))
            bond |= 1 << (L^2 + (x - 1) + L * (y - 1))
            push!(bonds, bond)
            # Add bond with the top-right spin
            bond = 0
            bond |= 1 << ((x - 1) + L * (y - 1))
            bond |= 1 << (L^2 + mod(x, L) + L * mod(y - 2, L))
            push!(bonds, bond)
            # Add bond with the bottom-right spin
            bond = 0
            bond |= 1 << ((x - 1) + L * (y - 1))
            bond |= 1 << (L^2 + mod(x, L) + L * (y - 1))
            push!(bonds, bond)
        end
    end
    return bonds
end

# The following function generates a list of all diagonal operators, with both integers.

# Given the size L, generate all diagonal terms.
function find_all_diagonals(L)
    plaquettes = find_all_plaquettes(L)
    vertices = find_all_vertices(L)
    bonds = find_all_bonds(L)
    diagonals = Vector{Int}[]
    for p in plaquettes
        push!(diagonals, [1, p])
    end
    for v in vertices
        push!(diagonals, [5, v])
    end
    for b in bonds
        push!(diagonals, [3, b])
    end
    return diagonals
end

# The following three functions are used for the QMC algorithm. The main function implements one 
# QMC sweep, consisting of all diagonal updates, off-diagonal updates, and spin flip updates. A 
# helper function is used to obtain the linked vertex representation of the configuration. Another 
# helper function is used to create an initial state where all the operators are padding identities.

# Given the current configuration, list of all diagonal operators, size L, inverse temperature beta,
# and the Hamiltonian parameters (J_1, J_2, C_1, C_2, C_3), perform one QMC sweep. 
function qmc_sweep(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, rng)
    # Diagonal updates
    M = length(operators)                   # Number of operators/states in the configuration
    N_d = length(diagonals)                 # Total number of diagonal operators available
    n = count(op -> op[1] != 0, operators)  # Number of non-padding operators in the configuration
    # Loop through the operators
    for i in 1:M
        op = operators[i]
        # Padding identity -> diagonal
        if op[1] == 0
            # Choose a random diagonal operator
            new_op = copy(rand(rng, diagonals))
            # Calculate the acceptance probability
            # ZZZZ (plaquette)
            if new_op[1] == 1
                p_acc = beta * N_d / (M - n) * (C_1 + (-1)^(count_ones(new_op[2] & states[i])))
            # ZZ/II (bond)
            elseif new_op[1] == 3
                if J_1 == 0
                    p_acc = beta * N_d / (M - n) * C_4
                else
                    p_acc = beta * N_d / (M - n) * J_1 * (C_2 + (-1)^(count_ones(new_op[2] & states[i])))
                end
            # IIII (vertex)
            else
                p_acc = beta * N_d / (M - n) * C_3
            end
            # Roll the dice
            if rand(rng, ) < p_acc
                # Update operator to diagonal operator
                operators[i] = new_op
                n += 1
            end
        # Diagonal -> padding identity
        elseif op[1] in (1, 3, 5)
            # Calculate acceptance probability
            # ZZZZ (plaquette)
            if op[1] == 1
                p_acc = (M - n + 1) / (beta * N_d * (C_1 + (-1)^(count_ones(op[2] & states[i]))))
            # ZZ/II (bond)
            elseif op[1] == 3
                if J_1 == 0
                    p_acc = (M - n + 1) / (beta * N_d * C_4)
                else
                    p_acc = (M - n + 1) / (beta * N_d * J_1 * (C_2 + (-1)^(count_ones(op[2] & states[i]))))
                end
            # IIII (vertex)
            else
                p_acc = (M - n + 1) / (beta * N_d * C_3)
            end
            # Roll the dice
            if rand(rng, ) < p_acc
                # Update operator to padding identity
                operators[i] = [0,0]
                n -= 1
            end
        end
    end

    # Off-diagonal updates
    # Find linked vertex representation
    # This has the form ((operator type, operator integer, operator index, [vertex index connected to leg 1, ...]), ...)
    linked_vertex = to_linked_vertex(states, operators, L)
    # Store a list of the operator sides to flip
    # Has the form: ((operator index, side), ...)
    vertices_to_flip = Set{Tuple{Int, Int}}()
    # Construct the clusters: We choose a leg as a starting point and traverse the list, keeping track
    # of visited legs. Then, we choose the next leg and repeat until all legs have been visited.
    # Has the form: ((vertex index, local leg number))
    visited_legs = Set{Tuple{Int, Int}}()
    for vertex_idx in 1:length(linked_vertex)
        for leg in 1:length(linked_vertex[vertex_idx][4])
            # Traverse all the vertices in the loop and calculate the acceptance probability
            stack = Tuple{Int, Int}[]                   # The stack has the form ((vertex index 1, leg number 1), ...)
            push!(stack, (vertex_idx, leg))
            cluster_vertices = Set{Tuple{Int, Int}}()   # This has the same form as vertices_to_flip but is specific to the cluster
            p_acc = 1
            while !isempty(stack)
                # Find the current vertex index and leg combination 
                curr_vertex_idx, curr_leg = pop!(stack)
                curr_vertex = linked_vertex[curr_vertex_idx]

                # Check that this has not been visited yet
                if (curr_vertex_idx, curr_leg) in visited_legs
                    continue
                end

                # Mark as visited
                push!(visited_legs, (curr_vertex_idx, curr_leg))

                # Determine if we are on the incoming or outgoing side (0 or 1)
                if curr_leg <= div(length(curr_vertex[4]), 2)
                    side = 0
                else
                    side = 1
                end

                # Add the neighbors to the stack
                # ZZZZ (plaquette) or ZZ (bond) if J_2 = 0
                if curr_vertex[1] == 1 || (curr_vertex[1] == 3 && J_2 == 0)
                    # Add the next leg, on both sides
                    if curr_leg < length(curr_vertex[4])
                        push!(stack, (curr_vertex_idx, curr_leg + 1))
                    else
                        push!(stack, (curr_vertex_idx, 1))
                    end
                    # Add the link
                    push!(stack, curr_vertex[4][curr_leg])
                # Bond or vertex
                else
                    # Add the next leg, same side
                    if side == 0
                        if curr_leg < div(length(curr_vertex[4]), 2)
                            push!(stack, (curr_vertex_idx, curr_leg + 1))
                        else
                            push!(stack, (curr_vertex_idx, 1))
                        end
                    else
                        if curr_leg < length(curr_vertex[4])
                            push!(stack, (curr_vertex_idx, curr_leg + 1))
                        else
                            push!(stack, (curr_vertex_idx, div(length(curr_vertex[4]), 2) + 1))
                        end
                    end
                    # Add the link
                    push!(stack, curr_vertex[4][curr_leg])
                end

                # Determine if the operator + side is new
                if !((curr_vertex[3], side) in cluster_vertices)
                    # Update acceptance probability
                    # First time encountering the vertex (across previous cluster and current cluster)
                    if !((curr_vertex[3], xor(side, 1)) in cluster_vertices) && !((curr_vertex[3], xor(side, 1)) in vertices_to_flip)
                        # ZZ/II -> XX
                        if curr_vertex[1] == 3 && J_2 != 0
                            # Evaluate the matrix element using the incoming state (either is fine)
                            state_idx = curr_vertex[3]
                            if J_1 == 0
                                p_acc /= C_4 / J_2
                            else
                                p_acc /= J_1 / J_2 * (C_2 + (-1)^(count_ones(curr_vertex[2] & states[state_idx])))
                            end
                        # XX -> ZZ/II
                        elseif curr_vertex[1] == 4
                            # Evaluate the matrix using the opposite (i.e. unchanged) of the specified side
                            state_idx = mod1((curr_vertex[3] + xor(side, 1)), M)
                            if J_1 == 0
                                p_acc *= C_4 / J_2
                            else
                                p_acc *= J_1 / J_2 * (C_2 + (-1)^(count_ones(curr_vertex[2] & states[state_idx])))
                            end
                        # IIII -> XXXX
                        elseif curr_vertex[1] == 5
                            p_acc /= C_3
                        # XXXX -> IIII
                        elseif curr_vertex[1] == 2
                            p_acc *= C_3
                        end
                    # Second time encountering the vertex, i.e. we need to reverse the flip
                    else
                        # ZZ/II -> XX -> ZZ/II
                        if curr_vertex[1] == 3 && J_2 != 0
                            # Evaluate the matrix element using the specified side, and wrap around
                            state_idx = mod1(curr_vertex[3] + side, M)
                            if J_1 == 0
                                p_acc *= C_4 / J_2
                            else
                                p_acc *= J_1 / J_2 * (C_2 + (-1)^(count_ones(curr_vertex[2] & states[state_idx])))
                            end
                        # XX -> ZZ/II -> XX
                        elseif curr_vertex[1] == 4
                            # Evaluate the matrix element using the specified side, and wrap around
                            state_idx = mod1(curr_vertex[3] + side, M)
                            if J_1 == 0
                                p_acc /= C_4 / J_2
                            else
                                p_acc /= J_1 / J_2 * (C_2 + (-1)^(count_ones(curr_vertex[2] & states[state_idx])))
                            end
                        # IIII -> XXXX -> IIII
                        elseif curr_vertex[1] == 5
                            p_acc *= C_3
                        # XXXX -> IIII -> XXXX
                        elseif curr_vertex[1] == 2
                            p_acc /= C_3
                        end
                    end
                end

                # Append the vertex number and side
                push!(cluster_vertices, (curr_vertex[3], side))
            end

            # Roll the dice
            if rand(rng, ) < 0.5         # Proposal probability
                if rand(rng, ) < p_acc   # Acceptance probability
                    union!(vertices_to_flip, cluster_vertices)
                end
            end
        end
    end

    # Flip the operators and reconstruct the spin states
    # Go from imaginary time = 1 to M, and flip the operators.
    to_flip = 0     # Keep track if we need to flip each spin
    found = 0       # Keep track if we have determined whether we need to flip each spin
    for i in 1:M
        # Flip the operator
        if xor(((i,0) in vertices_to_flip), ((i,1) in vertices_to_flip))
            # ZZ/II -> XX
            if operators[i][1] == 3
                operators[i][1] = 4
            # XX -> ZZ/II
            elseif operators[i][1] == 4
                operators[i][1] = 3
            # IIII -> XXXX
            elseif operators[i][1] == 5
                operators[i][1] = 2
            # XXXX -> IIII
            elseif operators[i][1] == 2
                operators[i][1] = 5
            end
        end
        # Update the list of first encounters
        if (i,0) in vertices_to_flip
            to_flip = xor(to_flip, (xor(((1 << (2*L^2)) - 1), found) & operators[i][2]))
            found = found | operators[i][2]
        end
        if (i,1) in vertices_to_flip
            found = found | operators[i][2]
        end
    end
    # Flip the initial state and reconstruct the rest
    states[1] = xor(states[1], to_flip)
    for i in 1:M-1
        # Flip if XX or XXXX, otherwise keep the same
        if operators[i][1] == 2 || operators[i][1] == 4
            states[i+1] = xor(operators[i][2], states[i])
        else
            states[i+1] = states[i]
        end
    end

    # Spin flip update
    # Loop through each spin
    for spin in 1:2*L^2
        # Loop through operators and update acceptance probability of flipping the whole column of spin
        p_acc = 1
        for i in 1:M
            # The weight only changes if there is a ZZ or ZZZZ operator acting on the spin
            op = operators[i]

            # Skip if padding identity
            if op[1] == 0
                continue
            end
            # Check if the operator acts on the spin
            if (op[2] & (1 << (spin - 1))) != 0
                # ZZZZ
                if op[1] == 1
                    p_acc *= (C_1 - (-1)^(count_ones(op[2] & states[i]))) / (C_1 + (-1)^(count_ones(op[2] & states[i])))
                # ZZ
                elseif op[1] == 3 && J_1 != 0
                    p_acc *= (C_2 - (-1)^(count_ones(op[2] & states[i]))) / (C_2 + (-1)^(count_ones(op[2] & states[i])))
                end
            end
        end
        # Roll the dice
        if rand(rng, ) < p_acc
            for i in 1:M
                states[i] = xor(states[i], (1 << (spin - 1)))
            end
        end
    end
    return states, operators
end

# Convert the list of states and operators to the linked vertex list. This is a list of vertices, with each vertex being a list containing
# four elements: the operator type, the operator integer, the operator index (imaginary time), and a list of the vertices (linked vertex 
# index) of the vertices attached to each leg. The convention is that leg 1 is the incoming spin with the smallest index, then we go through 
# the larger indices of the incoming spins, and repeat for the outgoing spins.
function to_linked_vertex(states, operators, L)
    # Create the vertices
    linked_vertex = Tuple{Int, Int, Int, Vector{Tuple{Int, Int}}}[]
    for (i, op) in enumerate(operators)
        # Determine if non-padding operator
        if op[1] != 0
            push!(linked_vertex, (op[1], op[2], i, [(0, 0) for _ in 1:(2 * count_ones(op[2]))]))
        end
    end
    # Determine the links by looping through the vertices and keeping track of which ones have most recently acted on each spin, as well as
    # the local leg number of the previous vertex. Keep track of the first instance to account for periodicity.
    most_recent = [(0, 0) for _ in 1:2*L^2]
    first_seen = [(0, 0) for _ in 1:2*L^2]
    for (v_idx, vertex) in enumerate(linked_vertex)
        # Determine the spins the vertex acts on
        spins = findall(i -> vertex[2] & (1 << (i - 1)) != 0, 1:ndigits(vertex[2], base=2))
        # For each spin, connect the incoming of the current vertex to the outgoing of the previous vertex
        for (leg, spin) in enumerate(spins)
            if most_recent[spin] != (0,0)
                # Connect current incoming to previous outgoing
                vertex[4][leg] = most_recent[spin]
                # Connect previous outgoing to current incoming
                linked_vertex[most_recent[spin][1]][4][most_recent[spin][2]] = (v_idx, leg)
            else
                first_seen[spin] = (v_idx, leg)
            end
            # Update the list
            most_recent[spin] = (v_idx, leg + length(spins))
        end
    end
    for spin in 1:2*L^2
        if first_seen[spin] != (0,0)
            # Connect incoming of first seen to outgoing of most recent
            linked_vertex[first_seen[spin][1]][4][first_seen[spin][2]] = most_recent[spin]
            # Connect outgoing of most recent to incoming of first seen
            linked_vertex[most_recent[spin][1]][4][most_recent[spin][2]] = first_seen[spin]
        end
    end
    return linked_vertex
end

# Given the configuration, measure all the necessary raw data
function measure_all(states, operators, L)
    # Measure number of operators (total and for each type)
    n = count(op -> op[1] != 0, operators)
    n_B = count(op -> op[1] == 1, operators)
    n_A = count(op -> op[1] == 2, operators)
    n_zz = count(op -> op[1] == 3, operators)
    n_xx = count(op -> op[1] == 4, operators)

    # Measure M_z, M_z^2 and M_z^4
    mz = 0
    mz_2 = 0
    mz_4 = 0
    for (state, op) in zip(states, operators)
        if op[1] != 0
            # Compute M_z
            n_down = count_ones(state)
            mz_p = 1 - 2 * n_down / (2*L^2)
            # Update
            mz += mz_p
            mz_2 += mz_p^2
            mz_4 += mz_p^4
        end
    end
    mz /= n
    mz_2 /= n
    mz_4 /= n

    # Measure structure factors
    S_0 = 0
    S_kmin = 0
    for (state, op) in zip(states, operators)
        if op[1] != 0
            # Compute structure factors for this state
            M_0 = 0
            M_kmin_x = 0
            M_kmin_y = 0
            # Loop through each spin
            for j in 0:2*L^2-1
                # Determine the x and y positions of the spin, with the origin at the top-left vertex
                # Horizontal link
                if j < L^2
                    r_x = mod(j, L) + 3/2
                    r_y = div(j, L) + 1
                # Vertical link
                else 
                    r_x = mod(j, L) + 1
                    r_y = div(j - L^2, L) + 3/2
                end
                # Calculate Z_j
                z_j = 1 - 2 * Int((state >> j) & 1)
                # Update the structure factors
                M_0 += z_j
                M_kmin_x += exp(im * 2*pi/L * r_x) * z_j
                M_kmin_y += exp(im * 2*pi/L * r_y) * z_j
            end
            S_0 += abs2(M_0)
            S_kmin += abs2(M_kmin_x) + abs2(M_kmin_y)
        end
    end
    S_0 /= n * (2*L^2)
    S_kmin /= 2*n * (2*L^2)

    return n, n_A, n_B, n_zz, n_xx, mz, mz_2, mz_4, S_0, S_kmin
end