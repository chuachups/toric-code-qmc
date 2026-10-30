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
# The cluster construction propogates the spin lines along all flippable operators (IIII, XXXX, ZZ/II (if J_2 =/= 0), XX)
# but passes straight through the non-flippable operators (ZZZZ, ZZ (if J_2 =/= 0)).
# There are three additional updates to help with ergodicity and convergence:
# 1. Spin column flips: Flip all spins at a certain site for all imaginary time
# 2. Dual loop updates: Flip all spins along a dual loop for all imaginary time
# 3. Interval updates: Similar to cluster update, but consists of only two flippable operators. The spin line passes through
#    all operators in between.

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
            # [REVISION 2026-10-23] The L^2 offset was on the wrong pair of spins: the top/bottom spins used the
            # horizontal-link formula and the left/right spins used the vertical-link formula. That made every
            # "vertex" equal to a plaquette for L=2, and made vertices anticommute with plaquettes for L>=3.
            # Corrected: top/bottom are vertical links (offset L^2), left/right are horizontal links (no offset).
            # Find the top spin (vertical link whose lower end is (x,y), i.e. the link above the vertex)
            vertex |= 1 << (L^2 + (x - 1) + L * mod(y - 2, L))
            # Find the bottom spin (vertical link below the vertex)
            vertex |= 1 << (L^2 + (x - 1) + L * (y - 1))
            # Find the left spin (horizontal link to the left of the vertex)
            vertex |= 1 << (mod(x - 2, L) + L * (y - 1))
            # Find the right spin (horizontal link to the right of the vertex)
            vertex |= 1 << ((x - 1) + L * (y - 1))
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

# [REVISION 2026-10-23] New helper for the topological-sector update (see REVISION_LOG.md).
# Given the size L, return a list of all non-contractible dual loops (X Wilson loops), each stored as an integer
# whose binary digits are 1 for the spins on the loop. A dual loop is a set of spins such that every plaquette
# contains an even number of them, so flipping the loop leaves every B_p (and every vertex term) unchanged, but
# it flips the conserved Z Wilson loop that it crosses. Candidates: the L horizontal links of one column, and the
# L vertical links of one row. Only candidates with even overlap with every plaquette are kept (this filter
# keeps the dual loops for either assignment of the plaquettes/vertices, e.g. in the X basis).
function find_all_dual_loops(L)
    plaquettes = find_all_plaquettes(L)
    candidates = Int[]
    for a in 1:L
        h_column = 0    # horizontal links (x = a, all y)
        h_row = 0       # horizontal links (all x, y = a)
        v_column = 0    # vertical links (x = a, all y)
        v_row = 0       # vertical links (all x, y = a)
        for b in 1:L
            h_column |= 1 << ((a - 1) + L * (b - 1))
            h_row |= 1 << ((b - 1) + L * (a - 1))
            v_column |= 1 << (L^2 + (a - 1) + L * (b - 1))
            v_row |= 1 << (L^2 + (b - 1) + L * (a - 1))
        end
        append!(candidates, [h_column, h_row, v_column, v_row])
    end
    loops = filter(c -> all(p -> iseven(count_ones(c & p)), plaquettes), unique(candidates))
    # There must be L loops in each of the two winding directions
    @assert length(loops) == 2 * L "Expected $(2L) dual loops, found $(length(loops))"
    return loops
end

# [REVISION 2026-10-23] New helper for the topological-sector update (see REVISION_LOG.md).
# Given the configuration and a dual loop, return the weight ratio W(loop flipped at all times) / W(current).
# Only diagonal operators with an odd overlap with the loop change their matrix element; for a true dual loop
# these are the ZZ bonds crossing it (ZZZZ always has even overlap, but is handled anyway for safety).
# XXXX, IIII, XX and II (J_1 = 0) operators have spin-independent matrix elements.
function dual_loop_ratio(states, operators, loop, J_1, C_1, C_2)
    ratio = 1.0
    for i in 1:length(operators)
        op = operators[i]
        if op[1] == 0 || iseven(count_ones(op[2] & loop))
            continue
        end
        # ZZZZ (plaquette)
        if op[1] == 1
            ratio *= (C_1 - (-1)^(count_ones(op[2] & states[i]))) / (C_1 + (-1)^(count_ones(op[2] & states[i])))
        # ZZ (bond); if J_1 = 0 this is II with constant weight C_4
        elseif op[1] == 3 && J_1 != 0
            ratio *= (C_2 - (-1)^(count_ones(op[2] & states[i]))) / (C_2 + (-1)^(count_ones(op[2] & states[i])))
        end
    end
    return ratio
end

# [REVISION 2026-10-24] Lattice lists are cached per L (they used to be rebuilt in every sweep).
const LATTICE_CACHE = Dict{Int, Any}()
function lattice_lists(L)
    get!(LATTICE_CACHE, L) do
        (vertices = find_all_vertices(L), bonds = find_all_bonds(L), dual_loops = find_all_dual_loops(L))
    end
end

# [REVISION 2026-10-23, optimized 2026-10-24] Local "interval" updates (star-interval and bond-interval), see
# REVISION_LOG.md Sec. 10-11.
# Why they are needed: in the cluster update a line that reaches a vertex operator (IIII/XXXX), or a bond operator
# when J_2 != 0, always takes all legs of that side together. A single spin line can never pass straight through
# such an operator, although that leaves its weight unchanged (its matrix element does not depend on the spins).
# Changes of that kind (e.g. creating a pair of XXXX on one star while vertex operators of the neighbouring stars
# sit in between) were therefore never proposed, which biased <A> (and made <XX> converge very slowly).
#
# The move: given a mask S (the 4 spins of a star, or the 2 spins of a bond) and an ordered pair of operators
# a != c that both act on exactly S, toggle both (IIII <-> XXXX, or ZZ/II <-> XX) and flip the spins of S on every
# time slice strictly after a and up to c (cyclically, so the interval may wrap around imaginary time).
# - Consistency: off-diagonal operators inside the interval map flipped states to flipped states; the toggled
#   endpoints absorb the flip (a now flips S on its way out, c on its way in).
# - Weight ratio: the two endpoint toggles, times every diagonal operator strictly inside the interval with an odd
#   overlap with S (ZZZZ: (C_1 - B)/(C_1 + B), ZZ: (C_2 - zz)/(C_2 + zz)). XXXX, IIII, XX and II have
#   spin-independent elements. An endpoint bond operator keeps its own Z_iZ_j (both of its spins flip).
# - Proposal ("pairs" option), both symmetric because toggling moves no operator and the same (a, c) undoes it:
#     :uniform (default) a and c independently uniform among the k operators on S (a == c: no move); the original
#              version. Long intervals are possible in one step, which decorrelates the slow XX / sector mode best.
#     :nearby  a uniform among the k operators on S, c = the next one on S after a (cyclically); short intervals,
#              higher acceptance, slightly cheaper, but in case 2 the XX error was about 2x larger at equal sweeps
#              (REVISION_LOG.md Sec. 11), so it is not the default.
#
# Cost (optimized version): one O(M) pass builds, for every mask, the list of its operators and, for every spin, the
# sorted list of operators acting on it, plus each operator's current eigenvalue sgn = B_p or Z_iZ_j. A proposal then
# walks only the operators on the spins of S inside the interval (each operator counted once), and an accepted move
# only flips the eigenvalues it walked. The per-slice states are rebuilt once at the end (states[1] is flipped by S
# for every accepted interval that wraps through slot 1), so the stage costs O(M) per sweep instead of
# O(number of masks * M) (review comment, 2026-10-24).
function interval_stage!(states, operators, L, J_1, J_2, C_1, C_2, C_3, C_4, rng; pairs = :uniform)
    M = length(operators)
    N = 2 * L^2
    lat = lattice_lists(L)
    do_bonds = J_2 != 0     # for J_2 = 0 there are no XX operators, so every bond toggle would have weight 0
    by_mask = Dict{Int, Vector{Int}}()
    by_spin = [Int[] for _ in 1:N]
    sgn = zeros(Int8, M)
    for i in 1:M
        op = operators[i]
        op[1] == 0 && continue
        if op[1] in (2, 5) || (do_bonds && op[1] in (3, 4))
            push!(get!(by_mask, op[2], Int[]), i)
        end
        if op[1] in (1, 3, 4)
            sgn[i] = iseven(count_ones(op[2] & states[i])) ? 1 : -1
        end
        m = op[2]
        while m != 0
            push!(by_spin[trailing_zeros(m) + 1], i)
            m &= m - 1
        end
    end
    # Weight of an endpoint operator of type t at position p (stars: IIII/XXXX; bonds: ZZ/II/XX)
    endpoint_weight(t, p) = t == 5 ? C_3 : t == 2 ? 1.0 : t == 4 ? J_2 : (J_1 == 0 ? C_4 : J_1 * (C_2 + sgn[p]))
    flip_slot1 = 0
    walked = Int[]
    masks = do_bonds ? vcat(lat.vertices, lat.bonds) : lat.vertices
    for S in masks
        idx = get(by_mask, S, nothing)
        (idx === nothing || length(idx) < 2) && continue
        if pairs == :nearby
            k = rand(rng, 1:length(idx))
            a, c = idx[k], idx[mod1(k + 1, length(idx))]
        else
            a, c = rand(rng, idx), rand(rng, idx)
            a == c && continue
        end
        ratio = 1.0
        for p in (a, c)
            t = operators[p][1]
            ratio *= endpoint_weight(toggled_type(t), p) / endpoint_weight(t, p)
        end
        # Operators strictly inside the cyclic interval (a, c), visited through the lists of the spins of S
        empty!(walked)
        m = S
        while m != 0
            j = trailing_zeros(m)
            v = by_spin[j + 1]
            ranges = a < c ? (searchsortedfirst(v, a + 1):(searchsortedfirst(v, c) - 1),) :
                             (searchsortedfirst(v, a + 1):length(v), 1:(searchsortedfirst(v, c) - 1))
            for r in ranges, q in r
                p = v[q]
                ov = operators[p][2] & S
                # count each operator once (at its lowest spin in S); only odd overlaps change anything
                (trailing_zeros(ov) == j && isodd(count_ones(ov))) || continue
                push!(walked, p)
                t = operators[p][1]
                if t == 1
                    ratio *= (C_1 - sgn[p]) / (C_1 + sgn[p])
                elseif t == 3 && J_1 != 0
                    ratio *= (C_2 - sgn[p]) / (C_2 + sgn[p])
                end
            end
            m &= m - 1
        end
        # Roll the dice
        if rand(rng, ) < ratio
            operators[a][1] = toggled_type(operators[a][1])
            operators[c][1] = toggled_type(operators[c][1])
            for p in walked
                sgn[p] = -sgn[p]        # eigenvalues of walked ZZZZ/ZZ/XX flip (sgn = 0 for the others)
            end
            if a > c
                flip_slot1 = xor(flip_slot1, S)
            end
        end
    end
    # Rebuild the per-slice states once
    states[1] = xor(states[1], flip_slot1)
    for i in 1:M-1
        if operators[i][1] == 2 || operators[i][1] == 4
            states[i+1] = xor(operators[i][2], states[i])
        else
            states[i+1] = states[i]
        end
    end
    return states, operators
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

# [REVISION 2026-10-23] New helpers for the exact cluster acceptance probability (see REVISION_LOG.md).
# Operator types are toggled by the cluster update as IIII <-> XXXX and ZZ/II <-> XX.
toggled_type(op_type) = op_type == 5 ? 2 : op_type == 2 ? 5 : op_type == 3 ? 4 : 3

# Matrix element of a "side-type" operator (IIII/XXXX, and ZZ/II/XX if J_2 != 0) of the given type. state is the
# original incoming state of the operator. Flipping a side of a bond operator flips both of its spins, so Z_iZ_j has
# the same value on both sides and in every flipped version, namely its value in the original incoming state.
function side_op_weight(op_type, op, state, J_1, J_2, C_2, C_3, C_4)
    if op_type == 5         # IIII
        return C_3
    elseif op_type == 2     # XXXX
        return 1.0
    elseif op_type == 4     # XX
        return J_2
    elseif J_1 == 0         # II
        return C_4
    else                    # ZZ
        return J_1 * (C_2 + (-1)^(count_ones(op[2] & state)))
    end
end

# Weight ratio W(after) / W(before) of the operator op at imaginary time t when the current cluster is flipped.
# "before" = original configuration + previously accepted clusters (vertices_to_flip); "after" = also the current
# cluster (cluster_vertices). state = states[t] is the ORIGINAL incoming state (states are updated only after all
# clusters have been decided).
# - Pass-through operators (ZZZZ, and ZZ/II if J_2 = 0) are stored with their leg numbers. Every spin line through
#   them flips one of their spins and contributes 2 legs, so the matrix element depends only on the parity of the
#   number of flipped lines. An II bond (J_1 = 0) has the spin-independent weight C_4.
# - Side-type operators are stored with sides 0 (incoming) and 1 (outgoing); the type is toggled iff exactly one side
#   is flipped.
# A zero new weight (possible for C_1 = 1 or C_2 = 1) gives a zero ratio instead of the old 0/0 = NaN.
function operator_flip_ratio(op, state, t, cluster_vertices, vertices_to_flip, J_1, J_2, C_1, C_2, C_3, C_4)
    if op[1] == 1 || (op[1] == 3 && J_2 == 0)
        if op[1] == 3 && J_1 == 0
            return 1.0
        end
        C = op[1] == 1 ? C_1 : C_2
        num_legs = 2 * count_ones(op[2])
        legs_before = count(j -> (t, j) in vertices_to_flip, 1:num_legs)
        legs_after = legs_before + count(j -> (t, j) in cluster_vertices, 1:num_legs)
        b = (-1)^(count_ones(op[2] & state))
        w_before = C + (iseven(div(legs_before, 2)) ? b : -b)
        w_after = C + (iseven(div(legs_after, 2)) ? b : -b)
        return w_after / w_before
    else
        f0_before = (t, 0) in vertices_to_flip
        f1_before = (t, 1) in vertices_to_flip
        f0_after = xor(f0_before, (t, 0) in cluster_vertices)
        f1_after = xor(f1_before, (t, 1) in cluster_vertices)
        type_before = xor(f0_before, f1_before) ? toggled_type(op[1]) : op[1]
        type_after = xor(f0_after, f1_after) ? toggled_type(op[1]) : op[1]
        return side_op_weight(type_after, op, state, J_1, J_2, C_2, C_3, C_4) /
               side_op_weight(type_before, op, state, J_1, J_2, C_2, C_3, C_4)
    end
end

# The following three functions are used for the QMC algorithm. The main function implements one 
# QMC sweep, consisting of all diagonal updates, off-diagonal updates, and spin flip updates. A 
# helper function is used to obtain the linked vertex representation of the configuration. Another 
# helper function is used to create an initial state where all the operators are padding identities.

# Given the current configuration, list of all diagonal operators, size L, inverse temperature beta,
# and the Hamiltonian parameters (J_1, J_2, C_1, C_2, C_3), perform one QMC sweep. 
# [REVISION 2026-10-24] New keyword interval_pairs (:uniform, default, or :nearby), the pair choice of the interval updates.
function qmc_sweep(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, rng, diagnostic; interval_pairs = :uniform)
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
                    # [REVISION 2026-10-23] II bond operators (weight C_4) must be removable for J_2 = 0 as well. The
                    # removal probability used to be 0 when J_1 = J_2 = 0, so II operators were inserted but never
                    # removed, and n (and M) grew without bound.
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
    # Find the linked vertex representation
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
            touched_ops = Set{Int}()                    # [REVISION 2026-10-23] imaginary times of all operators touched by the cluster
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
                    # Add the opposite leg
                    if curr_leg <= div(length((curr_vertex[4])), 2)
                        push!(stack, (curr_vertex_idx, curr_leg + div(length(curr_vertex[4]),2)))
                    else
                        push!(stack, (curr_vertex_idx, curr_leg - div(length(curr_vertex[4]),2)))
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

                # [REVISION 2026-10-23] The acceptance probability is no longer updated leg by leg here (the old running
                # multiply/divide gave 0/0 = NaN when C_1 = 1 or C_2 = 1, and used the ZZ formula for II bonds when
                # J_1 = J_2 = 0). Instead we record the operator, and compute its exact weight ratio once the cluster
                # is complete (see "Exact acceptance probability" below).
                push!(touched_ops, curr_vertex[3])

                # Append the vertex number and side/leg number
                if curr_vertex[1] == 1 || (curr_vertex[1] == 3 && J_2 == 0)
                    push!(cluster_vertices, (curr_vertex[3], curr_leg))
                else
                    push!(cluster_vertices, (curr_vertex[3], side))
                end
            end

            # [REVISION 2026-10-23] Exact acceptance probability of this cluster: the product over all touched operators
            # of W(operator after flipping this cluster) / W(operator now), where "now" already includes the flips of
            # previously accepted clusters (vertices_to_flip). The stored states are still the original ones.
            p_acc = 1.0
            for t in touched_ops
                p_acc *= operator_flip_ratio(operators[t], states[t], t, cluster_vertices, vertices_to_flip, J_1, J_2, C_1, C_2, C_3, C_4)
            end

            # Roll the dice
            if rand(rng, ) < 0.5         # Proposal probability
                if rand(rng, ) < p_acc   # Acceptance probability
                    union!(vertices_to_flip, cluster_vertices)
                end
            end

            # Diagnostics
            if diagnostic == true
                if length(cluster_vertices) != 0
                    println("Cluster with p_acc $p_acc")
                    for (t, label) in ((2, "XXXX"), (5, "IIII"), (4, "XX"), (3, J_1 == 0 ? "II" : "ZZ"))
                        unique_count = 0
                        duplicated_count = 0

                        for i in eachindex(operators)
                            operators[i][1] == t || continue
                            n = count(x -> x[1] == i, cluster_vertices)
                            unique_count += (n == 1)
                            duplicated_count += (n == 2)
                        end

                        println("$label: unique=$unique_count, duplicated=$duplicated_count")
                    end
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
            if operators[i][1] == 3 && J_2 != 0
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
        if !(operators[i][1] == 1 || (operators[i][1] == 3 && J_2 == 0))
            if (i,0) in vertices_to_flip
                to_flip = xor(to_flip, (xor(((1 << (2*L^2)) - 1), found) & operators[i][2]))
                found = found | operators[i][2]
            end
            if (i,1) in vertices_to_flip
                found = found | operators[i][2]
            end
        else
            if operators[i][1] == 1
                num_legs = 4
            else
                num_legs = 2
            end
            for j in 1:num_legs
                if (i,j) in vertices_to_flip
                    # Retain only the j-th 1 of the operator
                    spin_leg = operators[i][2]
                    for _ in 1:j-1
                        spin_leg &= spin_leg - 1
                    end
                    spin_leg = spin_leg & (-spin_leg)
                    to_flip = xor(to_flip, (xor(((1 << (2*L^2)) - 1), found) & spin_leg))
                    found = found | spin_leg
                end
            end
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

    # [REVISION 2026-10-23] Topological-sector (X Wilson loop) update, see REVISION_LOG.md.
    # The Z Wilson loops commute with every term when J_2 = 0, and the only other move that changes them (the
    # single-spin column flip above) has to create two fluxes along the whole imaginary-time axis, so at large
    # beta the sector never changed. Here we propose flipping all the spins of one non-contractible dual loop at
    # every imaginary time. All XXXX/XX operators still connect consecutive states, and every ZZZZ/IIII operator
    # keeps its matrix element (even overlap), so only operators with an odd overlap with the loop (the ZZ bonds
    # crossing it) change weight. Metropolis acceptance with the exact weight ratio.
    for loop in lattice_lists(L).dual_loops     # [REVISION 2026-10-24] cached instead of rebuilt every sweep
        p_acc = dual_loop_ratio(states, operators, loop, J_1, C_1, C_2)
        # Roll the dice
        if rand(rng, ) < p_acc
            for i in 1:M
                states[i] = xor(states[i], loop)
            end
        end
    end

    # [REVISION 2026-10-23, optimized 2026-10-24] Local interval updates (see the comment above interval_stage!):
    # one star-interval proposal per star and, for J_2 != 0, one bond-interval proposal per bond.
    # interval_pairs = :uniform (default) or :nearby.
    interval_stage!(states, operators, L, J_1, J_2, C_1, C_2, C_3, C_4, rng; pairs = interval_pairs)
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
    M = length(states)

    # Measure M_z, M_z^2 and M_z^4
    mz = 0
    mz_2 = 0
    mz_4 = 0
    for state in states
        # Compute M_z
        n_down = count_ones(state)
        mz_p = 1 - 2 * n_down / (2*L^2)
        # Update
        mz += mz_p
        mz_2 += mz_p^2
        mz_4 += mz_p^4
    end
    mz /= M
    mz_2 /= M
    mz_4 /= M

    # Measure structure factors
    S_kmin = 0
    for state in states
        # Compute structure factors for this state
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
            M_kmin_x += exp(im * 2*pi/L * r_x) * z_j
            M_kmin_y += exp(im * 2*pi/L * r_y) * z_j
        end
        S_kmin += abs2(M_kmin_x) + abs2(M_kmin_y)
    end
    S_kmin /= 2*M * (2*L^2)

    return n, n_A, n_B, n_zz, n_xx, mz, mz_2, mz_4, S_kmin
end

# [REVISION 2026-10-23] New measurement for the susceptibility estimator (see REVISION_LOG.md).
# Return the averages of M_z and M_z^2 (per qubit) over the n propagated states, i.e. the states just before each
# non-padding operator. The standard SSE estimator chi = beta*N*[(sum_p m_p)^2 + sum_p m_p^2] / (n(n+1)) is then
# beta*N*(n*m_n^2 + m_2_n)/(n+1). The averages of measure_all run over all M slices (padding included); that is fine
# for equal-time observables, but inserting them into this formula biased chi_M upwards (for large M, on average by
# beta*N*n*(m_2_n - m_n^2)/(n+1)^2). For n = 0 the single state is used (the formula then gives beta*N*m^2).
function measure_m_propagated(states, operators, L)
    N = 2 * L^2
    m_sum = 0.0
    m_2_sum = 0.0
    n = 0
    for i in eachindex(operators)
        if operators[i][1] != 0
            m_p = 1 - 2 * count_ones(states[i]) / N
            m_sum += m_p
            m_2_sum += m_p^2
            n += 1
        end
    end
    if n == 0
        m_p = 1 - 2 * count_ones(states[1]) / N
        return m_p, m_p^2
    end
    return m_sum / n, m_2_sum / n
end

# Given the current configuration, list of all diagonal operators, size L, inverse temperature beta,
# and the Hamiltonian parameters (J_1, J_2, C_1, C_2, C_3), perform one QMC sweep. 
# [REVISION 2026-10-23] Note: qmc_sweep_old is not called anywhere and was NOT revised (it still has the old bugs).
function qmc_sweep_old(states, operators, diagonals, L, beta, J_1, J_2, C_1, C_2, C_3, C_4, rng, diagnostic)
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

            # Diagnostics
            if diagnostic == true
                if length(cluster_vertices) != 0
                    println("Cluster with p_acc $p_acc")
                    for (t, label) in ((2, "XXXX"), (5, "IIII"), (4, "XX"), (3, J_1 == 0 ? "II" : "ZZ"))
                        unique_count = 0
                        duplicated_count = 0

                        for i in eachindex(operators)
                            operators[i][1] == t || continue
                            n = count(x -> x[1] == i, cluster_vertices)
                            unique_count += (n == 1)
                            duplicated_count += (n == 2)
                        end

                        println("$label: unique=$unique_count, duplicated=$duplicated_count")
                    end
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