# Summary of the revisions to the toric-code SSE package

This document is written for the author of the original code. It explains
what was changed in `toric-code-qmc - revised` relative to the original
package (`toric-code-qmc - Copy 2026-10-23`), and why. Each change is marked
in the source with a `[REVISION 2026-10-23]` or `[REVISION 2026-10-24]`
comment. `REVISION_LOG.md` is the detailed, chronological lab notebook (every
test and number). This file is the summary.

Scope of the review and of all tests: **L = 2**, mostly the Z basis, and the
two L = 2 cases flagged in `info.txt`:

- **case 1:** β = 5, J₁ = 0.1, J₂ = 0
- **case 2:** β = 5, J₁ = J₂ = 0.1

---

## 1. The short version

| # | Problem in the original code | Effect | Change |
|---|---|---|---|
| 1 | **Wrong star (vertex) operators** in `find_all_vertices` (Julia and Python) | Not the toric code. At L = 2 every "star" equals a plaquette; at L ≥ 3 stars anticommute with plaquettes. ED used the same masks, so ED vs SSE could not detect it. | Fixed masks. Stale cache files are detected and regenerated. |
| 2 | **Topological sectors frozen** (Z basis) | Each run stays in the sector of its random initial state. This is the main cause of the wrong values in case 1. | New X-Wilson-loop (dual-loop) update. |
| 3 | **Cluster update cannot pass a single spin line through a vertex (or bond) operator** | ⟨A⟩ (and E, ZZ, M², …) systematically biased in every seed (case 1); XX in case 2 converges extremely slowly. | New star-interval and bond-interval updates. |
| 4 | Smaller bugs | Thermalization crash, unbounded n for J₁ = J₂ = 0, NaN for C = 1, biased χ_M, import side effect, … | Fixed (Sec. 5). |

With all revisions, cases 1 and 2 agree with exact diagonalization (Sec. 7).

**Important for reading the old results:** every individual update of the
original code (diagonal update, cluster update, column flip) satisfied
detailed balance. This was checked by brute force. The wrong numbers came
from the model geometry (problem 1) and from the chain being unable to reach
parts of configuration space (problems 2 and 3), not from wrong acceptance
probabilities. The per-run jackknife error bars could not detect this,
because each run is self-consistent in the region it samples. Comparing many
independent seeds reveals it immediately.

---

## 2. Problem 1: the star operators (model geometry)

**Code convention** (Julia, 1-based x, y):

- horizontal link `i = (x-1) + L(y-1)` (vertex (x,y) to its left);
- vertical link `i = L² + (x-1) + L(y-1)` (vertex (x,y) above it).

In the original `find_all_vertices`, the L² offset was on the wrong pair of
spins. The "top/bottom" spins used the horizontal formula, and "left/right"
used the vertical formula. At L = 2 this gives

```
"stars"    : [0,2,4,5] [0,2,6,7] [1,3,4,5] [1,3,6,7]
plaquettes : [0,2,4,5] [0,2,6,7] [1,3,4,5] [1,3,6,7]     <- identical
```

So the simulated Hamiltonian was −Σ_S X_S − Σ_S Z_S − … with X and Z on the
*same* four sets. At L ≥ 3 these "stars" anticommute with plaquettes (18
pairs at L = 3), so the model is not even a commuting stabilizer model there.

**Why the ED comparison never showed it:** `exact_diagonalization.py`
contains the same formula, so ED and SSE simulated the same (wrong) model. At
L = 2 the full-trace thermal averages of the wrong model also happen to be
equal (E, A, B, ZZ, M², …) or nearly equal (χ_z in case 2: 1.613205 vs
1.613185) to those of the correct one at the tested points. The sector
structure, however, is different.

**Fix** (`Julia/SSE.jl` and `Python/exact_diagonalization.py`):

```
top    = L^2 + (x-1) + L*mod(y-2,L)   # vertical link above the vertex
bottom = L^2 + (x-1) + L*(y-1)        # vertical link below
left   = mod(x-2,L) + L*(y-1)         # horizontal link to the left
right  = (x-1) + L*(y-1)              # horizontal link to the right
```

After the fix, at L = 2 the stars are [0,1,4,6], [2,3,4,6], [0,1,5,7],
[2,3,5,7]. At L = 2, 3, 4: all stars commute with all plaquettes, and every
ZZ bond lies inside exactly one star and anticommutes with exactly two. The
self-duality needed by the X-basis mode (a qubit relabelling that maps stars
↔ plaquettes and bonds → bonds) still holds; this was checked numerically.

**Cache** (`Julia/run_job.jl`): `Precomputed/diagonals_L*.h5` stored the old
masks. The diagonal-operator list is now always built from the lattice
functions. A cached file that differs is overwritten with a warning, so an old
cache cannot silently bring the bug back.

---

## 3. Problem 2: frozen topological sectors → X-Wilson-loop update

**Why sectors freeze.** On the torus, a product of Z along a non-contractible
loop (a Z Wilson loop; at L = 2 for example Z₀Z₁ and Z₄Z₆) commutes with every
star, plaquette and ZZ term. For J₂ = 0 it therefore commutes with H, and the
trace splits into 2² = 4 topological sectors. ED averages over all four: in
case 1 their weights are 0.609, 0.130, 0.130, 0.130. In the Z-basis SSE,
XXXX and all diagonal operators conserve these loops. The only original move
that changes them is the single-spin column flip, which creates two fluxes
along the whole imaginary-time axis. Its acceptance is about
(1/3)^(number of ZZZZ operators on that spin), roughly 10⁻¹⁴ at β = 5. In
tests: **0 sector changes in 2·10⁵ sweeps.** Case 1 therefore converged to
the values of whichever sector the random initial state was in.

**The new move** (end of `qmc_sweep`; helpers `find_all_dual_loops`,
`dual_loop_ratio`):

1. For each of the 2L non-contractible **dual** loops ℓ, propose flipping
   the spins of ℓ on **every** time slice. A dual loop is a set of spins that
   every plaquette contains an even number of times; at L = 2 they are
   {0,2}, {1,3}, {4,5}, {6,7}.
2. Accept with probability min(1, R), where R is the product, over ZZ
   operators whose bond crosses the loop, of (C₂ − zz)/(C₂ + zz).

Why this works:

- X_ℓ commutes with every B_p and every star, so the configuration stays
  valid, and ZZZZ, IIII and XXXX keep their matrix elements.
- X_ℓ anticommutes with the Z Wilson loop it crosses, so the sector changes
  **without creating fluxes**. The acceptance is O(1) at small J₁.
- The proposal does not depend on the configuration and undoes itself, so
  Metropolis acceptance is exact.

This is not an operator insertion. The operator string is untouched; it is
equivalent to changing the trace state |α(0)⟩ → X_ℓ|α(0)⟩. In case 1 the
sector now changes about every 4–5 sweeps (mean acceptance about 9% per
proposal, 4 proposals per sweep).

---

## 4. Problem 3: missing local moves → star-interval and bond-interval updates

**What was missing.** Your cluster update treats ZZZZ (and ZZ when J₂ = 0) as
*pass-through*: a spin line enters and leaves on the same spin. It treats
IIII/XXXX (and ZZ/XX when J₂ ≠ 0) as *junctions*: a line that reaches one
always pulls in all legs on that side. Every cluster flip, and every
combination of them, therefore flips each side of a junction operator either
completely or not at all.

The SSE weight allows more than that. The matrix elements of IIII, XXXX and XX
do not depend on the spins. So flipping the incoming **and** outgoing leg of
one spin of such an operator (a line passing straight through it) is also a
valid change with unchanged weight. The natural local move, "create a pair of
XXXX on star s at τ₁ < τ₂ and flip the 4 spins of s in between", crosses the
vertex operators of neighbouring stars (which share 2 spins with s) in exactly
this way. Your cluster update can never propose it, and the diagonal update
cannot create or remove XXXX at all.

Note that XXXX operators *were* created and removed constantly: vertex-toggling
clusters are accepted about 84% of the time. The problem is which
*arrangements* can be reached. The result was a seed-independent bias:
⟨A⟩ = 0.917 vs exact 0.956 in case 1, with E, ZZ, M² and the sector weights
off by 15–23σ. In case 2 the same mechanism for bond operators made XX
converge extremely slowly.

**The new moves** (`interval_stage!` in `SSE.jl`, called at the end of
`qmc_sweep`):

1. Take a mask S: the 4 spins of a star (star-interval), or the 2 spins of a
   NN bond (bond-interval, only for J₂ ≠ 0).
2. Pick an ordered pair a ≠ c of operators acting on exactly S (IIII/XXXX for
   stars, ZZ/II/XX for bonds).
3. Toggle both (IIII ↔ XXXX, ZZ/II ↔ XX), and flip the spins of S on every
   time slice strictly after a and up to c (cyclically).
4. Weight ratio: the two endpoint toggles, times every diagonal operator
   strictly inside the interval with an *odd* overlap with S. For ZZZZ the
   factor is (C₁ − B)/(C₁ + B); for ZZ it is (C₂ − zz)/(C₂ + zz). For a star
   this includes only ZZ bonds crossing the star boundary. For a bond it
   includes the two plaquettes that contain one of its spins, so the move
   creates a pair of fluxes for that interval, exactly as a physical XX pair
   does.
5. The proposal is symmetric, because toggling moves no operator and the same
   (a, c) undoes the move. Metropolis acceptance is therefore exact.

**Pair choice:** `qmc_sweep(...; interval_pairs = :uniform)` is the default:
a and c are independently uniform among the operators on S.
`interval_pairs = :nearby` takes c as the next operator on S after a. It gives
shorter intervals and slightly lower cost, but in case 2 it gave about twice
the XX error at equal sweeps, so it is not the default.

**Cost:** one O(M) pass builds per-mask and per-spin operator lists and each
operator's current eigenvalue. A proposal then walks only the operators on the
spins of S, and `states` is rebuilt once at the end. The whole stage costs a
constant 8–10% of a sweep, at all L tested (2–5).

---

## 5. Smaller bug fixes

| File | Problem | Fix |
|---|---|---|
| `Julia/thermalize.jl` | New cutoff slots were padded with `[0]` instead of `[0, 0]`, so `reduce(hcat, operators)` in `write_output_file` threw `DimensionMismatch`. The new slots also copied `states[end]` (the state *before* the last operator), which is wrong when that operator is off-diagonal. | Pad with `[0, 0]`; the new slots get `states[1]` (the state after the last operator, by periodicity). |
| `Julia/SSE.jl` (diagonal update) | For J₁ = J₂ = 0 the removal probability of II bond operators was 0, so they accumulated and n and M grew without bound (pure toric-code runs were wrong and slow). | II operators can always be removed with probability (M − n + 1)/(β N_d C₄). |
| `Julia/SSE.jl` (cluster acceptance) | The acceptance was built leg by leg with a running multiply and divide. This (a) gives 0/0 = NaN when C₁ = 1 or C₂ = 1, and (b) applied the ZZ formula to II bonds at J₁ = J₂ = 0, whose weight C₄ does not depend on the spins. | After the cluster is complete, the exact ratio W(after)/W(before) is computed for each touched operator (`operator_flip_ratio`, `side_op_weight`, `toggled_type`). For C = 2 this is numerically identical to the old value: the same seeds give bit-identical runs. |
| `Julia/SSE.jl`, `thermalize.jl`, `simulate.jl`, `write_file.jl`, `Python/data_analysis.py` | χ_M was computed from averages over **all M slices** (padding included) inserted into the n-operator formula. That is not the standard estimator; it is biased upwards (small here). | New `measure_m_propagated` gives ⟨m⟩ and ⟨m²⟩ over the n propagated states. They are stored as `measurements/m_n` and `measurements/m_2_n`, and `data_analysis.py` uses them: χ = βN[(Σ m_p)² + Σ m_p²]/(n(n+1)). Old files without these datasets fall back to the old estimator, with a warning. |
| `Python/data_analysis.py` | `print_all(...)` ran at import time, so `from data_analysis import ...` (as in `sse_vs_ed.py`) crashed unless that particular file existed. | Wrapped in `if __name__ == "__main__":`. |
| `Julia/SSE.jl` (`qmc_sweep_old`) | Unused. | A comment says it was **not** revised (it still contains the old bugs). |

---

## 6. Reference: what is new or changed, file by file

**`Julia/SSE.jl`**

- `find_all_vertices`: fixed (Sec. 2).
- New `find_all_dual_loops`, `dual_loop_ratio`, `lattice_lists` (cached
  stars, bonds and dual loops per L).
- New `interval_stage!` (star- and bond-interval updates, Sec. 4).
- New `toggled_type`, `side_op_weight`, `operator_flip_ratio` (exact cluster
  acceptance, Sec. 5).
- `qmc_sweep`: the II removal fix; the exact cluster acceptance; then, after
  the original column flips, the dual-loop update and the interval stage.
  New keyword `interval_pairs = :uniform`. **The positional signature is
  unchanged.**
- New `measure_m_propagated` (for χ_M).
- `to_linked_vertex` and `measure_all`: unchanged.

**Other files**

- `Julia/thermalize.jl`: padding fix; records `m_n`, `m_2_n`.
- `Julia/simulate.jl`: records `m_n`, `m_2_n`.
- `Julia/write_file.jl`: new keyword arguments `m_n_arr`, `m_2_n_arr`
  (written only if non-empty). The positional signature is unchanged.
- `Julia/run_job.jl`: cache check for `Precomputed/`.
- `Python/exact_diagonalization.py`: `find_all_vertices` fixed.
- `Python/data_analysis.py`: χ_M estimator; import guard.
- **Unchanged:** `main.jl`, `make_jobs.py`, `jobs.py`, `run_local.py`,
  `run_pbs.py`, `sse_vs_ed.py`, `thermalization_analysis.py`. The constants
  chosen by `make_jobs.py` (C₁ = C₂ = 2, C₃ = 1, C₄) are unchanged.

**Compatibility.** Old HDF5 files can still be analysed (with the χ_M
warning). Results of the old and new code are not comparable, because the
model geometry changed (except for the coincidental L = 2 full-trace
equality).

---

## 7. How the revised code was validated

1. **Brute-force detailed balance.** For every move (cluster, dual loop,
   star interval, bond interval), the code's acceptance ratio was compared
   with W(after)/W(before) computed from the full configuration, and the
   flip application with an independent implementation. Tens of thousands of
   proposals at many parameter points, including C = 1, J₁ = 0, J₂ ≠ 0:
   0 mismatches.
2. **Exact test of the complete Markov chain.** With a small fixed cutoff M,
   the SSE samples exactly the truncated series Z_M = Σ_{n≤M} βⁿ/n! Tr Tⁿ
   (T = C − H), whose expectations follow from matrix powers. This test is
   independent of how fast the chain mixes. All moves together reproduce ⟨n⟩,
   the operator counts, the sector probabilities and ⟨m²⟩ to about 10⁻⁴, at
   three parameter points with J₂ ≠ 0, for both pair options.
3. **Full production path.** `main.jl` thermalize → simulate (warm start) →
   HDF5 → `data_analysis.py`, in the Z and X bases.
4. **Cases vs ED.** ED was recomputed with the corrected masks and checked
   against an independent dense ED. Mean ± standard error over independent
   seeds, 10⁵ sweeps each:

**Case 1 (β = 5, J₁ = 0.1, J₂ = 0), 8 seeds**

| | ED | original code | revised |
|---|---|---|---|
| E | −8.1876 | −8.158 ± 0.039 | −8.1881 ± 0.0036 |
| A | 0.9565 | 0.895 ± 0.011 | 0.9562 ± 0.0006 |
| ZZ | 0.2261 | 0.361 ± 0.052 | 0.2263 ± 0.0003 |
| M² | 0.3188 | 0.434 ± 0.062 | 0.3196 ± 0.0005 |
| sector weights | 0.609 / 0.130 / 0.130 / 0.130 | frozen (0 or 1 per run) | 0.609 / 0.131 / 0.130 / 0.130 |

**Case 2 (β = 5, J₁ = J₂ = 0.1), 28 seeds**

| | ED | revised |
|---|---|---|
| E | −8.3379 | −8.3357 ± 0.0039 |
| A / B | 0.9623 / 0.9623 | 0.9623 ± 0.0004 / 0.9633 ± 0.0009 |
| ZZ / XX | 0.2000 / 0.2000 | 0.2016 ± 0.0008 / 0.1939 ± 0.0042 |
| M² | 0.2937 | 0.2952 ± 0.0007 |
| χ_M | 1.6132 | 1.6222 ± 0.0060 |

In the original code, case 2 did not settle down at all: XX wandered from
0.05 to 0.87 and M² from 0.01 to 0.52 between blocks of 5 000 sweeps.

---

## 8. Known limitations and advice

- **System size:** spin configurations are one `Int64` bitmask, so the code
  works only for N = 2L² ≤ 63, i.e. **L ≤ 5** (it fails at L = 6). Larger L
  needs a multi-word representation.
- **L = 3 (third case of `info.txt`):** tested after the revision
  (`REVISION_LOG.md` Sec. 13.1). At β = 1, J₁ = J₂ = 1, 16 seeds × 10⁶ sweeps
  agree with exact diagonalization within 0.9σ for every observable. The chain
  mixes slowly there: A, B, ZZ, XX and M² have τ_int ≈ 600–1000 sweeps, from a
  slow exchange of weight between the X-type and Z-type terms. That is why
  short runs look unconverged and the M² trace shows long plateaus. Use
  ≥ 10⁶ sweeps per seed at this point.
- **Not tested:** the X basis away from J₁ = J₂. At J₁ = J₂ the X-basis run
  is the same chain as the Z-basis run, so the X-basis bookkeeping (J₁ ↔ J₂
  swap, A ↔ B, ZZ ↔ XX, C₄) has not been checked against ED for J₁ ≠ J₂.
- **Case 2:** the ~2σ tension of the 28-seed test was a statistical
  fluctuation (`REVISION_LOG.md` Sec. 13.2). With 32 seeds × 10⁶ sweeps,
  every observable and P(Z-loop 00) agree with ED within 0.9σ, and the other
  three sector weights (16 seeds) within 1.6σ. The slowest mode at this point is the XX count
  (τ_int ≈ 460 sweeps). The Z-loop sector itself decorrelates in ≈ 15 sweeps.
- **Error bars:** always run several independent seeds (≥ 8) from different
  random initial states, and compare the spread between seeds with the
  in-run error bars. Plot block time series of M², XX and the Z-loop sector
  (for J₂ = 0 this also shows immediately whether the sector moves).
