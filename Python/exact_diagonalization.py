# The Quantum Monte Carlo (QMC) is based on a Stochasitic Series Expansion (SSE) of the partition function
# in the computational basis. Each term in the expansion consists of a basis state and a list of terms from
# the Hamiltonian (vertex, plaquette, bond, or identity), which may repeat. The following functions are used
# to enumerate all the vertex, plaquette, and bond terms. This is also useful for Exact Diagonalization (ED).

# The spins on the lattice can be indexed 0, 1, 2, ..., 2L^2 - 1.
# The following convention will be used:
# - For spins on horizontal links: i(x, y) = x + Ly where (x,y) is the vertex to the left of the link.
# - For spins on vertical links: i(x, y) = L^2 + x + Ly where (x,y) is the vertex above the link.

import numpy as np

# This function returns a list of all the vertex operators, each stored as an integer whose binary digits
# are 1 for the spin indices
def find_all_vertices(L):
    vertices = []
    # Loop through all the vertices (x,y)
    for x in range(L):
        for y in range(L):
            vertex = 0
            # Find the top spin
            vertex |= 1 << (x + L * ((y - 1) % L))
            # Find the bottom spin
            vertex |= 1 << (x + L * y)
            # Find the left spin
            vertex |= 1 << (L**2 + (x - 1) % L + L * y)
            # Find the right spin
            vertex |= 1 << (L**2 + x + L * y)
            vertices.append(vertex)
    return vertices

# This function returns a list of all the plaquette operators, each stored as an integer whose binary digits
# are 1 for the spin indices
def find_all_plaquettes(L):
    plaquettes = []
    # Loop through all the plaquettes (x,y)
    for x in range(L):
        for y in range(L):
            plaquette = 0
            # Find the top spin
            plaquette |= 1 << (x + L * y)
            # Find the bottom spin
            plaquette |= 1 << (x + L * ((y + 1) % L))
            # Find the left spin
            plaquette |= 1 << (L**2 + x + L * y)
            # Find the right spin
            plaquette |= 1 << (L**2 + (x + 1) % L + L * y)
            plaquettes.append(plaquette)
    return plaquettes

# This function returns a list of all the nearest neighbor bonds, each stored as an integer whose binary digits
# are 1 for the spin indices
def find_all_bonds(L):
    bonds = []
    # Loop through all the horizontal links
    for x in range(L):
        for y in range(L):
            # Add bond with top-left spin
            bond = 0
            bond |= 1 << (x + L * y)
            bond |= 1 << (L**2 + x + L * ((y - 1) % L))
            bonds.append(bond)
            # Add bond with the bottom-left spin
            bond = 0
            bond |= 1 << (x + L * y)
            bond |= 1 << (L**2 + x + L * y)
            bonds.append(bond)
            # Add bond with the top-right spin
            bond = 0
            bond |= 1 << (x + L * y)
            bond |= 1 << (L**2 + (x + 1) % L + L * ((y - 1) % L))
            bonds.append(bond)
            # Add bond with the bottom-right spin
            bond = 0
            bond |= 1 << (x + L * y)
            bond |= 1 << (L**2 + (x + 1) % L + L * y)
            bonds.append(bond)
    return bonds

# Exact Diagonalization
# The Hamiltonian commutes with P_Z (= Z_1 Z_2 ... Z_2L^2) and P_X (= X_1 X_2 ... X_2L^2)
# The Hamiltonian also commutes with T_x (which shifts all spins horizontally) and T_y (vertically)
# The goal is to block diagonaloze the Hamiltonian first, then use numpy to diagonalize each block

# Construct the block diagonalized matrix
# The logic is as follows:
# - Each spin has an index from 0 to 2L^2 - 2. We can assign a binary number to each computational basis state,
#   with 0 (1) in bit n representing a spin +1 (-1) for the spin with index n. Then, the basis states can be
#   numbered from 0 to 2^(2L^2) - 1. Note that these basis states are already eigenstates of P_Z.
# - Afterwards, we diagonalize the states according to P_X. We can take each state and add (or subtract)
#   its "partner" (|a> +- P_x |a>). To avoid double counting, we can loop through only the first half of the 
#   computational basis states.
# - Next, we diagonalize the states according to T_x. For each p_x eigenstate |a>, we can form a momentum-x
#   eigenstate with eigenvalue k_x by repeatedly adding exp(-i k_x r) T_x^r |a> for r = 0, 1, ..., L - 1.
#   Some of these eigenstates may be zero, or non-unique. We keep only the unique ones.
# - Finally, we can do the same but for T_y, starting with each momentum-x eigenstate.
# - Now, we have different sectors of eigenstates with different p_z, k_x, k_y, and p_x. Since the corresponding
#   operators commute with the Hamiltonian, each sector block-diagonalizes the Hamiltonian. We simply have to
#   evaluate the matrix elements using each sector as a basis.

# To save on runtime, we pre-compute all the sectors at the start (for fixed L) so that we do not have to redo it
# for each new Hamiltonian (different values of J_1, J_2).
def find_sectors(L):
    # Each state will be stored as dictionary with the compuational basis states (integers) as keys, and the 
    # coefficients (complex numbers) as the values. The eigenstates in a sector are stored in a list, and all the
    # sector lists are stored in a dictionary with the keys being a tuple of the eigenvalues (p_z, p_x, k_x, k_y)
    # Intermediate dictionaries are used throughout, labeled by the latest eigenvalue which has been found.

    # Diagonalize P_X using the compuational basis states
    px_eigenvals_states = {}
    for i in range(2**(2*L**2-1)):
        # Determine the p_z eigenvalue
        p_z = 1 if i.bit_count() % 2 == 0 else -1

        # Loop through the p_x eigenvalues
        partner = (~i) & (2**(2*L**2) - 1)
        for p_x in [1, -1]:
            eigenstate = {i: 1/np.sqrt(2), partner: p_x/np.sqrt(2)}
            if (p_z, p_x) in px_eigenvals_states:
                px_eigenvals_states[(p_z, p_x)].append(eigenstate)
            else:
                px_eigenvals_states[(p_z, p_x)] = [eigenstate]

    # Diagonalize T_x using the p_x states
    kx_eigenvals_states = {}
    for n in range(L):
        # Loop through k_x = n * 2*pi/L where n = 0, 1, ..., L-1
        k_x = 2 * np.pi * n / L
        
        for eigenvalues, px_eigenstates in px_eigenvals_states.items():
            p_z = eigenvalues[0]
            p_x = eigenvalues[1]

            seen = set()   # Store the computational basis states that we have seen
            for px_eigenstate in px_eigenstates:
                # Construct the eigenstate with definite k_x (and p_z, p_x)
                eigenstate = diagonalize_k(px_eigenstate, k_x, L, "x")

                # Compute normalization constant
                N = 0
                for coef in eigenstate.values():
                    N += coef.real**2 + coef.imag**2

                # If zero, skip
                if N > 0:
                    # Normalize
                    N = np.sqrt(N)
                    for cbasis_state in eigenstate:
                        eigenstate[cbasis_state] /= N

                    # Determine if the state is new, by checking if there are any shared computational basis states
                    # zero coefficient (entire state is automatically identical by construction) among states in the
                    # same sector
                    if (p_z, p_x, k_x) in kx_eigenvals_states:
                        to_check = next(iter(eigenstate))
                        if to_check not in seen:
                            kx_eigenvals_states[(p_z, p_x, k_x)].append(eigenstate)
                            seen.update(eigenstate)
                    else:
                        kx_eigenvals_states[(p_z, p_x, k_x)] = [eigenstate]
                        seen.update(eigenstate)

    # Diagonalize T_y using the k_x eigenstates
    sectors = {}
    for n in range(L):
        # Loop through k_y = n * 2*pi/L where n = 0, 1, ..., L-1
        k_y = 2 * np.pi * n / L

        for eigenvalues, kx_eigenstates in kx_eigenvals_states.items():
            p_z = eigenvalues[0]
            p_x = eigenvalues[1]
            k_x = eigenvalues[2]

            seen = set() # Store the computational basis states that we have seen
            for kx_eigenstate in kx_eigenstates:
                # Construct the eigenstate with definite k_y (and p_z, k_x)
                eigenstate = diagonalize_k(kx_eigenstate, k_y, L, "y")

                # Compute normalization constant
                N = 0
                for coef in eigenstate.values():
                    N += coef.real**2 + coef.imag**2

                # If zero, skip
                if N > 0:
                    # Normalize
                    N = np.sqrt(N)
                    for cbasis_state in eigenstate:
                        eigenstate[cbasis_state] /= N

                    # Determine if the state is new, by checking if there are any shared computational basis states
                    # zero coefficient (entire state is automatically identical by construction) among states in the
                    # same sector
                    if (p_z, p_x, k_x, k_y) in sectors:
                        to_check = next(iter(eigenstate))
                        if to_check not in seen:
                            sectors[(p_z, p_x, k_x, k_y)].append(eigenstate)
                            seen.update(eigenstate)
                    else:
                        sectors[(p_z, p_x, k_x, k_y)] = [eigenstate]
                        seen.update(eigenstate)

    print("Done computing symmetry sectors.")
    return sectors

# Now that all the basis states have been found, it is time to find the matrix elements.
# The full Hamiltonian is given by H = -H_tc - J_1 H_zz - J_2 H_xx where H_tc = H_a + H_b is the
# toric code Hamiltonian and H_zz and H_xx are the Ising perturbations with Z and X bonds respectively.
# Each of these 4 terms will be block diagonalized separately, and then later combined for custom
# values of J_1 an J_2. This needs to be done only once for each L, saving runtime.
# The magnetization operators are also in block form, and precomputed here.
def block_H(L):
    # Pre-compute sectors, vertices, plaquettes, and bonds
    sectors = find_sectors(L)
    vertices = find_all_vertices(L)
    plaquettes = find_all_plaquettes(L)
    bonds = find_all_bonds(L)

    # Store the blocks as a list lists of numpy arrays
    op_blks = [] # Store H_a, H_b, H_zz, mz_2, mz_4, mx_2, mx_4
    mz_blks = {} # Store mz
    mx_blks = {} # Store mx

    # For the structure factors, find the fourier transform phase of each link position
    rx = np.zeros(2*L**2)
    ry = np.zeros(2*L**2)
    for j in range(2*L**2):
        if j < L**2:
            rx[j] = j % L + 3/2
            ry[j] = j // L + 1
        else:
            rx[j] = j % L + 1
            ry[j] = (j - L**2) // L + 3/2
    phase_x = np.exp(1j * 2 * np.pi / L * rx)   # Due to kmin along x
    phase_y = np.exp(1j * 2 * np.pi / L * ry)   # Due to kmin along y

    cbasis_to_eigenstate_dict = {}
    # Loop through each sector of eigenstates
    for sector, eigenstates in sectors.items():
        H_a = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex)
        H_b = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex)
        H_zz = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex) 
        H_xx = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex)
        mz_2 = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex)
        mx_2 = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex)
        S_kmin_z = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex)
        S_kmin_x = np.zeros((len(eigenstates), len(eigenstates)), dtype=complex)

        # (For later) Determine which eigenstate which computational basis component belongs to
        cbasis_to_eigenstate = {}
        for i, eigenstate in enumerate(eigenstates):
            for cbasis_state in eigenstate:
                cbasis_to_eigenstate[cbasis_state] = i
        cbasis_to_eigenstate_dict[sector] = cbasis_to_eigenstate # We will reuse this for mz and mx

        # Loop through the eigenstates
        for j in range(len(eigenstates)):
            ket = eigenstates[j]

            # We want to evaluate (c_1*<i_0| + ...)(H)(c_1|j_0> + ...)
            # To do this, we first act H on the ket. The resulting ket is a
            # linear combination of other eigenstates in the same sector, 
            # we simply have to find out which ones. 
            # By construction, the computational basis states are disjoint 
            # among the eigenstates (within a sector). We take advantage of this
            # to find out which eigenstates are present and their coefficients

            # Note: For the diagonal operators, the symmetry basis states have the same
            # matrix element as any one of the computational basis states which make it up
            # (This also requires that the operators commute with the symmetry operators.)
            cbasis_state = next(iter(ket))
            # Hamiltonian terms
            H_b[j, j] = sum(find_Z(plaquette, cbasis_state) for plaquette in plaquettes)
            H_zz[j, j] = sum(find_Z(bond, cbasis_state) for bond in bonds)
            # Magnetization squared (also k=0 structure factor)
            mz_2[j, j] = (1 - cbasis_state.bit_count() / L**2)**2
            # kmin structure factor
            M_kx = 0j
            M_ky = 0j
            for spin in range(2*L**2):
                # Find the spin of this basis state
                z_j = 1 - 2 * ((cbasis_state >> spin) & 1)
                # Add
                M_kx += phase_x[spin] * z_j
                M_ky += phase_y[spin] * z_j
            S_kmin_z[j, j] = (abs(M_kx)**2 + abs(M_ky)**2) / (4 * L**2)

            # For off diagonal operators
            # Determine the action on each computational basis state
            for cbasis_state in ket:
                # Act each off-diagonal term on the basis state
                for vertex in vertices:
                    cbasis_new = apply_X(vertex, cbasis_state)
                    # This term only contributes to <i|H|j> where |i> is the unique eigenstate that contains
                    # |cbasis_new>. This can be found using the dictionary from earlier. Similarly, we perform
                    # an inner product <i|j_0> where |j_0> is the computational basis state being considered
                    # then multiply by the coefficients. 
                    # Note: It is possible that X|j_0> belongs to a different sector. Since the overall H preserves
                    # which sector the eigenstate belongs to, we know the contributions will cancel out and can
                    # just assign 0.
                    if cbasis_new in cbasis_to_eigenstate:
                        i = cbasis_to_eigenstate[cbasis_new]
                        H_a[i, j] += np.conjugate(eigenstates[i][cbasis_new]) * ket[cbasis_state]
                for bond in bonds:
                    cbasis_new = apply_X(bond, cbasis_state)
                    # Same logic as the above.
                    if cbasis_new in cbasis_to_eigenstate:
                        i = cbasis_to_eigenstate[cbasis_new]
                        H_xx[i, j] += np.conjugate(eigenstates[i][cbasis_new]) * ket[cbasis_state]
                for spin in range(2*L**2):
                    for spin_2 in range(2*L**2):
                        cbasis_new = apply_X((1 << spin) ^ (1 << spin_2), cbasis_state)
                        # Same logic as above
                        if cbasis_new in cbasis_to_eigenstate:
                            i = cbasis_to_eigenstate[cbasis_new]
                            mx_2[i, j] += np.conjugate(eigenstates[i][cbasis_new]) * ket[cbasis_state] / (2*L**2)**2
                            S_kmin_x[i, j] += np.conjugate(eigenstates[i][cbasis_new]) * ket[cbasis_state] / (2*L**2) * (np.conjugate(phase_x[spin])) * phase_x[spin_2] / 2  # kmin along x
                            S_kmin_x[i, j] += np.conjugate(eigenstates[i][cbasis_new]) * ket[cbasis_state] / (2*L**2) * (np.conjugate(phase_y[spin])) * phase_y[spin_2] / 2  # kmin along y
                            

        op_blks.append([H_a, H_b, H_zz, H_xx, mz_2, mz_2 @ mz_2, mx_2, mx_2 @ mx_2, S_kmin_z, S_kmin_x])

    # Loop through each sector again to fill mz and mx
    # Note that these operators are not block diagonal, but are still in block form
    # mz
    seen_z = set()
    for sector, eigenstates in sectors.items():
        if sector in seen_z:
            continue

        # We must calculate the cross block <sector partner|mz (mx)|sector> where the partner sector has the opposite P_x (P_z)
        sector_partner_z = (sector[0], -sector[1] , sector[2], sector[3])
        cbasis_to_eigenstate_z = cbasis_to_eigenstate_dict[sector_partner_z]
        eigenstates_partner_z = sectors[sector_partner_z]

        seen_z.update([sector, sector_partner_z])

        # Repeat the same logic except we use the dictionary from the partner sector to determine the partner element
        mz = np.zeros((len(eigenstates_partner_z), len(eigenstates)), dtype=complex)

        # Loop through each symmetry eigenstate
        for j in range(len(eigenstates)):
            ket = eigenstates[j]
            # Determine the action on each element
            for cbasis_state in ket:
                if cbasis_state in cbasis_to_eigenstate_z:
                    i = cbasis_to_eigenstate_z[cbasis_state]
                    mz[i, j] += np.conjugate(eigenstates_partner_z[i][cbasis_state]) * (1 - cbasis_state.bit_count() / L**2) * ket[cbasis_state]
        
        mz_blks[(sector_partner_z, sector)] = mz
    # mx
    seen_x = set()
    for sector, eigenstates in sectors.items():
        if sector in seen_x:
            continue

        # We must calculate the cross block <sector partner|mz (mx)|sector> where the partner sector has the opposite P_x (P_z)
        sector_partner_x = (-sector[0], sector[1] , sector[2], sector[3])
        cbasis_to_eigenstate_x = cbasis_to_eigenstate_dict[sector_partner_x]
        eigenstates_partner_x = sectors[sector_partner_x]

        seen_x.update([sector, sector_partner_x])

        # Repeat the same logic except we use the dictionary from the partner sector to determine the partner element
        mx = np.zeros((len(eigenstates_partner_x), len(eigenstates)), dtype=complex)

        # Loop through each symmetry eigenstate
        for j in range(len(eigenstates)):
            ket = eigenstates[j]
            # Determine the action on each element
            for cbasis_state in ket:
                for spin in range(2*L**2):
                    cbasis_new = apply_X(1 << spin, cbasis_state)
                    # Same logic as above
                    if cbasis_new in cbasis_to_eigenstate_x:
                        i = cbasis_to_eigenstate_x[cbasis_new]
                        mx[i, j] += np.conjugate(eigenstates_partner_x[i][cbasis_new]) * ket[cbasis_state] / (2*L**2)

        mx_blks[(sector_partner_x, sector)] = mx

    return op_blks, mz_blks, mx_blks, sectors   

# Given a computational basis state (integer), apply X operators to a list of spins
def apply_X(X_string, state):
    return X_string ^ state # bitwise XOR flips the spins

# Given a computational basis state (integer), apply a translation operator T (x or y)
def apply_T(state, L, direction):
    translated = 0

    # Loop through each link (integer representation)
    for i in range(2*L**2):
        # Find the x, y from the integer representation
        offset = (i // L**2) * L**2
        local_i = i - offset
        x = local_i % L
        y = local_i // L
        
        # Determine new index of the spin after translation
        if direction == "x":
            x_new = (x + 1) % L
            y_new = y
        elif direction == "y":
            x_new = x
            y_new = (y + 1) % L
        i_new = offset + x_new + L * y_new

        # Assign the current value to the translated position of the new state
        if (state >> i) & 1:
            translated |= 1 << i_new

    return translated

# Given a computational basis state (integer), determint the eigenvalue of a product of Z operators
def find_Z(Z_string, state):
    return 1 if (Z_string & state).bit_count() % 2 == 0 else -1 # bitwise AND ensures only the 1's from the desired spins are counted

# Given a reference state (dictionary of computational basis state : coefficient pairs), k eigenvalue, and direction
# (x or y), find an eigenstate through repeated application of the translation operator
def diagonalize_k(ref_state, k, L, direction):
    translated_state = ref_state
    eigenstate = {}
    for r in range(L):
        # Add the current state to the eigenstate by looping through each computational basis state
        # and adding its coefficient * exp(-ikr)
        phase = np.exp(-1j * k * r)
        for cbasis_state in translated_state:
            if cbasis_state in eigenstate:
                eigenstate[cbasis_state] += translated_state[cbasis_state] * phase
            else:
                eigenstate[cbasis_state] = translated_state[cbasis_state] * phase
        # Translate the state
        translated_state_new = {}
        for cbasis_state in translated_state:
            translated_state_new[apply_T(cbasis_state, L, direction)] = translated_state[cbasis_state]
        translated_state = translated_state_new

    # Remove any computational basis states with zero coefficient
    for cbasis_state in list(eigenstate):
        if abs(eigenstate[cbasis_state]) < 1e-12:
            del eigenstate[cbasis_state]
    return eigenstate

# The goal is to use ED as a benchmark for SSE. We will compare the following observables at
# different values of (J_1, J_2, beta):
# - Energy <E>, Heat capacity C
# - <A>, <B>
# - <ZZ>, <XX>
# - <M_z^2>, <M_z^4> and chi_z
# - <M_x^2>, <M_x^4> and chi_x
# - chi_J

# This function finds all the observables
def find_all_observables(op_blks, mz_blks, mx_blks, sectors, L, J_1, J_2, beta):
    # Initialize the observables
    E = 0
    E_2 = 0
    A = 0
    B = 0
    ZZ = 0
    XX = 0
    mz_2 = 0
    mz_4 = 0
    mx_2 = 0
    mx_4 = 0
    S_kmin_z = 0
    S_kmin_x = 0
    V = 0
    chi_z = 0
    chi_x = 0
    chi_j = 0
    Z = 0

    # 1. Diagonalize and find all the energies, eigenstates
    # 2. Rewrite the relevant observables (blocks) in the energy eigenbasis
    # 3. Update observables
    r = J_2 / J_1
    spectrum = {}           # Dictionary of all the energies per sector
    eigenstate_blocks = {}  # Dictionary of matrices, each matrix contains the eigenstates in each sector, in the symmetry basis
    for op_blk, sector, in zip(op_blks, sectors.keys()):
        # Find energy eigenbasis
        H_a, H_b, H_zz, H_xx, mz_2_mat, mz_4_mat, mx_2_mat, mx_4_mat, S_kmin_z_mat, S_kmin_x_mat = op_blk
        H = -H_a - H_b - J_1*H_zz - J_2*H_xx
        E_blk, U_blk = np.linalg.eigh(H)
        spectrum[sector] = E_blk
        eigenstate_blocks[sector] = U_blk

        # Rewrite the relevant observables in the energy eigenbasis
        A_blk = U_blk.T.conj() @ H_a @ U_blk / L**2
        B_blk = U_blk.T.conj() @ H_b @ U_blk / L**2
        ZZ_blk = U_blk.T.conj() @ H_zz @ U_blk / (4*L**2)
        XX_blk = U_blk.T.conj() @ H_xx @ U_blk / (4*L**2)
        mz_2_blk = U_blk.T.conj() @ mz_2_mat @ U_blk
        mz_4_blk = U_blk.T.conj() @ mz_4_mat @ U_blk
        mx_2_blk = U_blk.T.conj() @ mx_2_mat @ U_blk
        mx_4_blk = U_blk.T.conj() @ mx_4_mat @ U_blk
        S_kmin_z_blk = U_blk.T.conj() @ S_kmin_z_mat @ U_blk
        S_kmin_x_blk = U_blk.T.conj() @ S_kmin_x_mat @ U_blk
        V_blk = 4*L**2 * (ZZ_blk + r * XX_blk)

        for i in range(len(E_blk)):
            E_i = E_blk[i]
            E += E_i * np.exp(-beta * E_i)
            E_2 += E_i**2 * np.exp(-beta * E_i)
            A += A_blk[i,i].real * np.exp(-beta * E_i)
            B += B_blk[i,i].real * np.exp(-beta * E_i)
            ZZ += ZZ_blk[i,i].real * np.exp(-beta * E_i)
            XX += XX_blk[i,i].real * np.exp(-beta * E_i)
            mz_2 += mz_2_blk[i, i].real * np.exp(-beta * E_i)
            mz_4 += mz_4_blk[i, i].real * np.exp(-beta * E_i)
            mx_2 += mx_2_blk[i, i].real * np.exp(-beta * E_i)
            mx_4 += mx_4_blk[i, i].real * np.exp(-beta * E_i)
            S_kmin_z += S_kmin_z_blk[i, i].real * np.exp(-beta * E_i)
            S_kmin_x += S_kmin_x_blk[i, i].real * np.exp(-beta * E_i)
            V += V_blk[i, i].real * np.exp(-beta * E_i)
            Z += np.exp(-beta * E_i)

            for j in range(len(E_blk)):
                E_j = E_blk[j]
                if E_i != E_j:
                    chi_j += (np.exp(-beta * E_j) - np.exp(-beta * E_i)) / (E_i - E_j) * abs(V_blk[i,j])**2
                else:
                    chi_j += beta * np.exp(-beta * E_i) * abs(V_blk[i,j])**2

    # Express mz and mx in the energy eigenbasis as well, and compute the magnetic susceptibilities
    for sector1, sector2 in mz_blks.keys():
        mz_blk = eigenstate_blocks[sector1].T.conj() @ mz_blks[(sector1, sector2)] @ eigenstate_blocks[sector2]
        # Update susceptibility
        E_blk1 = spectrum[sector1]
        E_blk2 = spectrum[sector2]
        for i in range(len(E_blk1)):
            E_i = E_blk1[i]
            for j in range(len(E_blk2)):
                E_j = E_blk2[j]
                if abs(E_i - E_j) > 10**(-12):
                    chi_z += (np.exp(-beta * E_j) - np.exp(-beta * E_i)) / (E_i - E_j) * abs(mz_blk[i,j])**2
                    chi_z += (np.exp(-beta * E_i) - np.exp(-beta * E_j)) / (E_j - E_i) * abs(mz_blk[i,j])**2
                else:
                    chi_z += beta * np.exp(-beta * E_i) * abs(mz_blk[i,j])**2
                    chi_z += beta * np.exp(-beta * E_j) * abs(mz_blk[i,j])**2
    for sector1, sector2 in mx_blks.keys():
        mx_blk = eigenstate_blocks[sector1].T.conj() @ mx_blks[(sector1, sector2)] @ eigenstate_blocks[sector2]
        # Update susceptibility
        E_blk1 = spectrum[sector1]
        E_blk2 = spectrum[sector2]
        for i in range(len(E_blk1)):
            E_i = E_blk1[i]
            for j in range(len(E_blk2)):
                E_j = E_blk2[j]
                if abs(E_i - E_j) > 10**(-12):
                    chi_x += (np.exp(-beta * E_j) - np.exp(-beta * E_i)) / (E_i - E_j) * abs(mx_blk[i,j])**2
                    chi_x += (np.exp(-beta * E_i) - np.exp(-beta * E_j)) / (E_j - E_i) * abs(mx_blk[i,j])**2
                else:
                    chi_x += beta * np.exp(-beta * E_i) * abs(mx_blk[i,j])**2
                    chi_x += beta * np.exp(-beta * E_j) * abs(mx_blk[i,j])**2

    # Normalize
    E /= Z
    E_2 /= Z
    A /= Z
    B /= Z
    ZZ /= Z
    XX /= Z
    mz_2 /= Z
    mz_4 /= Z
    mx_2 /= Z
    mx_4 /= Z
    S_kmin_z /= Z
    S_0_z = 2*L**2 * mz_2
    if S_0_z > S_kmin_z:
        xi_z = 1 / (2 * np.sin(np.pi / L)) * np.sqrt(S_0_z / S_kmin_z - 1)
    else:
        xi_z = 0
    S_kmin_x /= Z
    S_0_x = 2*L**2 * mx_2
    if S_0_x > S_kmin_x:
        xi_x = 1 / (2 * np.sin(np.pi / L)) * np.sqrt(S_0_x / S_kmin_x - 1)
    else:
        xi_x = 0
    V /= Z
    C = beta**2 * (E_2 - E**2)
    chi_z *= 2*L**2 / Z
    chi_x *= 2*L**2 / Z
    chi_j /= 2*L**2 * Z
    chi_j -= beta * V**2 / (2*L**2)

    # Pack into dictionary and return
    ed = {
        "E": E,
        "mz2": mz_2,
        "mz4": mz_4,
        "mx2": mx_2,
        "mx4": mx_4,
        "C": C,
        "chi_z": chi_z,
        "chi_x": chi_x,
        "chi_j": chi_j,
        "xi_z": xi_z,
        "xi_x": xi_x,
        "A": A,
        "B": B,
        "ZZ": ZZ,
        "XX": XX,
    }
    return ed