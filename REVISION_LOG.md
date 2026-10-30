# Revision log: `toric-code-qmc - revised`

This folder started as a copy of `toric-code-qmc - Copy 2026-10-23` (2026-10-23).
Every code change is marked in the source with a `[REVISION 2026-10-23]`
comment. The background analysis is in `../QMC_REVIEW_FINDINGS.md`. The test
scripts are in `../review_scripts/revised/`.

## 1. Summary of revisions

| # | File | Change |
|---|---|---|
| 1 | `Julia/SSE.jl` → `find_all_vertices` | **Star operators fixed.** The L² (vertical-link) offset was on the left/right spins instead of the top/bottom spins. Before the fix, every "star" was a plaquette at L = 2, and stars anticommuted with plaquettes at L ≥ 3. |
| 2 | `Python/exact_diagonalization.py` → `find_all_vertices` | Same fix, so that ED and SSE describe the same (now correct) model. |
| 3 | `Precomputed/diagonals_L*.h5` | **Deleted.** These cached the old vertex masks; `run_job.jl` regenerates them on first use. |
| 4 | `Julia/SSE.jl` → new `find_all_dual_loops(L)` | Lists the 2L non-contractible dual loops (X Wilson loops). |
| 5 | `Julia/SSE.jl` → new `dual_loop_ratio(...)` | Exact weight ratio for flipping one dual loop at all imaginary times. |
| 6 | `Julia/SSE.jl` → end of `qmc_sweep` | **New topological-sector update**: one Metropolis proposal per dual loop per sweep. |

**Second round of revisions (small fixes, details in Sec. 7):**

| # | File | Change |
|---|---|---|
| 7 | `Julia/thermalize.jl` (cutoff update) | Padding operators are now `[0, 0]` (was `[0]`, which crashed `write_output_file`). New slots get the state *after* the last operator, `states[1]` (was `states[end]`, the state *before* it). |
| 8 | `Julia/SSE.jl` (diagonal update) | II bond operators (J₁ = 0) can be removed for J₂ = 0 too. The removal probability was 0, so n and M grew without bound. |
| 9 | `Julia/SSE.jl` (cluster update) + new `toggled_type`, `side_op_weight`, `operator_flip_ratio` | The cluster acceptance probability is now computed exactly, once per cluster, as a product of W(after)/W(before) over the touched operators. This replaces the leg-by-leg multiply/divide, which (a) gave 0/0 = NaN for C₁ = 1 or C₂ = 1, and (b) applied the ZZ weight formula to II bonds (J₁ = J₂ = 0), whose weight C₄ does not depend on the spins. |
| 10 | `Julia/SSE.jl`, new `measure_m_propagated` | ⟨m⟩ and ⟨m²⟩ averaged over the n propagated states, for χ_M. |
| 11 | `Julia/thermalize.jl`, `Julia/simulate.jl`, `Julia/write_file.jl` | Record and store these as `measurements/m_n` and `measurements/m_2_n` (new keyword arguments of `write_output_file`). |
| 12 | `Python/data_analysis.py` | χ_M uses `m_n`/`m_2_n`, i.e. the standard estimator βN[(Σ m_p)² + Σ m_p²]/(n(n+1)). The old M-slice averages are used only as a fallback for old files, with a warning. The module-level `print_all(...)` call now runs only under `if __name__ == "__main__":`, so importing the module (as `sse_vs_ed.py` does) no longer crashes. |
| 13 | `Julia/run_job.jl` | The diagonal-operator list is always rebuilt from the lattice functions. A cached `Precomputed/diagonals_L*.h5` that differs (e.g. one written by the original code with the wrong star masks) is overwritten with a warning instead of being used silently. |
| 14 | `Julia/SSE.jl` (`qmc_sweep_old`) | A comment only: this function is unused and was **not** revised. |

No other files were changed.

Check of the corrected geometry, at L = 2, 3, 4: 0 star–plaquette
anticommutations; every ZZ bond anticommutes with exactly 2 stars and lies
inside exactly 1 star. At L = 2 the stars are now
[0,1,4,6], [2,3,4,6], [0,1,5,7], [2,3,5,7] and the plaquettes are
[0,2,4,5], [0,2,6,7], [1,3,4,5], [1,3,6,7].

Note: at L = 2 the old (wrong) model happens to have the same full-trace
thermal averages as the corrected one at the tested point (E, A, ZZ, M², M⁴ are
identical). That is why the ED-vs-SSE comparison never exposed the geometry
bug. The sector structure is different, though, and L ≥ 3 is not equivalent at
all.

## 2. Strategy of the topological-sector update

**Why it is needed.** On the torus, the Z Wilson loops (a product of Z along a
non-contractible loop of the primal lattice; at L = 2 for example Z₀Z₁ and
Z₄Z₆) commute with every star, plaquette and ZZ term. For J₂ = 0 they also
commute with the whole Hamiltonian. The trace therefore splits into 4
topological sectors (2^2 for the two winding directions), and ED averages over
all of them with their Boltzmann weights. In the Z-basis SSE, the vertex
(XXXX) operators and every diagonal operator conserve these loops. The only
original move that changes them is the single-spin column flip. That flip
changes two plaquettes for the whole imaginary-time axis, so its acceptance is
about (1/3)^(number of ZZZZ operators on the spin), roughly 10⁻¹⁴ at β = 5.
In practice the sector was frozen at its random initial value (0 changes in
2 × 10⁵ sweeps).

**The move.** Pick one non-contractible *dual* loop ℓ: a set of spins that
every plaquette contains an even number of times. On the code's lattice these
are the L horizontal links of one column, or the L vertical links of one row
(`find_all_dual_loops` builds these candidates and keeps those with even
overlap with all plaquettes, so it also works if the roles of the
plaquette/vertex lists are swapped). Propose X_ℓ on **every** time slice:
`states[i] ⊻= ℓ` for all i.

- **Consistency:** every off-diagonal operator (XXXX, XX) still maps
  `states[i]` to `states[i+1]`, because XOR with ℓ commutes with XOR with the
  operator mask. Diagonal operators stay diagonal.
- **Weight ratio:** only diagonal operators with an *odd* overlap with ℓ change
  their matrix element. For a dual loop, ZZZZ always has even overlap, so the
  only affected operators are the ZZ bonds that cross the loop. Each
  contributes (C₂ − zz)/(C₂ + zz), where zz is its current value. XXXX, IIII,
  XX and II have spin-independent elements.
  `dual_loop_ratio` also handles ZZZZ in case of odd overlap, for safety.
- **Effect:** X_ℓ anticommutes with the Z Wilson loop that it crosses once,
  so the move flips that loop's eigenvalue and changes the topological sector.
  B_p and the vertex terms are untouched, so no fluxes are created and the
  acceptance stays O(1) at small J₁.
- **Detailed balance:** the proposal is deterministic and self-inverse (the
  same ℓ again undoes it), and it does not depend on the configuration.
  Metropolis acceptance min(1, ratio) therefore satisfies detailed balance.
- **Cost:** O(M) per loop, 2L loops per sweep. This is negligible next to the
  cluster update.

## 3. Verification (scripts in `../review_scripts/revised/`)

| Test | Result |
|---|---|
| `test_loop_move.jl`: `dual_loop_ratio` vs brute-force W(flipped)/W(current) over the full configuration | 0 mismatches in 24 000 proposals: (β, J₁, J₂) = (5, 0.1, 0), (5, 0.1, 0.1), (1, 1, 0) |
| `test_apply.jl`, `test_apply_j2.jl`: unchanged cluster update re-checked with the new geometry | 0 mismatches (3 858 + 5 273 clusters); flip application identical in 600 sweeps |
| Fixed Python ED vs independent dense ED (`dense_ed.py`) | identical to all printed digits |

## 4. Test of case 1 (L = 2, β = 5, J₁ = 0.1, J₂ = 0)

ED of the corrected model, by sector (`true_sector_ed.py`): sector
(Z₀Z₁, Z₄Z₆) = (0,0) has weight 0.6087; (0,1), (1,0) and (1,1) have 0.1304
each.

All SSE runs use 8 independent seeds (random initial state), 20 000
thermalization + 100 000 measurement sweeps, and fixed M = 300 (⟨n⟩ ≈ 117).
Errors are the standard error over seeds.

| Observable | ED | Original code (sector frozen) | Revised (stars + loop move) | Revised + star-interval move (experimental, Sec. 5) |
|---|---|---|---|---|
| E | −8.18757 | −8.158 ± 0.039 | −8.0985 ± 0.0041 (+22σ) | **−8.1881 ± 0.0036** (−0.1σ) |
| C | 0.5655 | — | 1.23 ± 0.16 | **0.53 ± 0.17** |
| A | 0.95646 | 0.895 ± 0.011 | 0.9170 ± 0.0017 (−23σ) | **0.95619 ± 0.00055** (−0.5σ) |
| B | 1.0 | 1.0005 ± 0.0007 | 1.0006 ± 0.0004 | **1.0002 ± 0.0011** |
| ZZ | 0.22608 | 0.361 ± 0.052 | 0.2666 ± 0.0024 (+17σ) | **0.22630 ± 0.00031** (+0.7σ) |
| M² | 0.31882 | 0.434 ± 0.062 | 0.3503 ± 0.0019 (+17σ) | **0.31964 ± 0.00053** (+1.5σ) |
| M⁴ | 0.29437 | — | 0.3283 ± 0.0021 (+17σ) | **0.29517 ± 0.00056** (+1.4σ) |
| P(sector 00) | 0.6087 | 0 or 1 per run | 0.6488 ± 0.0027 (+15σ) | **0.6092 ± 0.0010** (+0.5σ) |
| P(sector 01) / (10) / (11) | 0.1304 each | 0 or 1 per run | 0.118 / 0.117 / 0.116 | **0.1311 / 0.1300 / 0.1298** |
| sector changes per 10⁵ sweeps | — | 0 | ≈ 21 600 | ≈ 23 600 |

(The "original" column uses the old masks, whose sector labels differ; see
`../QMC_REVIEW_FINDINGS.md` Sec. 4.)

**Reading the table:**

1. The loop move does what it should. Sectors now change about every 5 sweeps,
   and the seed-to-seed scatter that came from the random initial sector is
   gone.
2. A systematic bias remains with only revisions 1–6 (A, E, ZZ, M² and the
   sector weights are off by 15–23σ, identically in every seed). This is
   mechanism R2 of the findings file. With the experimental star-interval move
   added, **every observable agrees with ED within 1.5σ.** That confirms R2 is
   a missing-move (ergodicity) problem in the original update scheme, not a
   wrong acceptance formula.

## 5. Remaining issue and proposed next revision (not applied to the package yet)

**R2: why the original moves do not sample the XXXX operators correctly.**
(Corrected 2026-10-23. An earlier version of this paragraph claimed that
vertex-toggling clusters are almost always rejected. Measurement shows the
opposite. With `measure_moves.jl` at case 1: all 19 904 clusters toggle at
least one IIII/XXXX, typically 2, and their mean acceptance is 0.84. XXXX
operators are created and removed constantly.)

The real restriction is the cluster rule at vertex operators. In `qmc_sweep`,
a line that reaches an IIII/XXXX operator always pulls in **all 4 legs on that
side** ("next leg, same side"). Only ZZ/ZZZZ let a single line pass straight
through. Every cluster flip, and every product of cluster flips, therefore
flips each vertex-operator side either completely or not at all. But the SSE
weight allows more. Because the matrix element of IIII or XXXX does not
depend on the spins, flipping the incoming **and** outgoing leg of the same
spin (a line passing straight through) is also a valid change with an
unchanged weight. Example: a star-interval flip of star s between its
operators a and b crosses vertex operators of the neighbouring stars that
share 2 spins with s. At each such operator c, the line passes straight
through on 2 of c's 4 legs, which is a "mixed" side. No combination of cluster
flips can produce that, so such configuration changes are never proposed. The
diagonal update cannot help directly either, because it can only remove IIII,
never XXXX.

The result is a systematic, seed-independent bias: A is too low by about 4% at
case 1, and the bias is larger at large β and small J₁. I haven't determined
whether this is strict non-ergodicity or only an extremely slow route; the
tight agreement between seeds (A = 0.908–0.922) suggests the former.

**Fix that was tested:** `../review_scripts/revised/star_interval_addon.jl`
(`star_interval_ratio`, `star_interval_apply!`, `star_interval_update!`). For
each star s, pick an ordered pair (a, b) of vertex operators acting on s.
Toggle both (IIII ↔ XXXX), and flip the 4 spins of s on every time slice
strictly after a and up to b (cyclically). XXXX operators in between still
connect consistent states. ZZZZ has even overlap with a star, so only ZZ bonds
inside the interval with odd overlap with s change weight. The proposal is
self-inverse and symmetric (the number of vertex operators at s is unchanged),
so Metropolis acceptance is exact. This move supplies exactly the missing
"line passes straight through a vertex operator" changes. A more systematic
alternative is to give vertex operators a stochastic two-way rule in the
cluster construction itself (junction vs straight-through, chosen by heat
bath), similar to bond choices in Swendsen–Wang. Brute-force check: 0 mismatches in 17 610
proposals (`test_addon.jl`). It was called once per star per sweep after
`qmc_sweep` (`case1_rev_addon.jl`).

**To integrate it into the package:** move the three functions into `SSE.jl`
and call `star_interval_update!` at the end of `qmc_sweep`. For J₂ ≠ 0 the
same construction works. XX operators in the interval have spin-independent
weights and stay consistent, but bond operators are then side-type. This needs
its own brute-force test before production use.

## 6. Status of the issues listed after the first revision

| Issue | Status |
|---|---|
| `thermalize.jl`: `[0]` padding and the wrong appended state | **fixed** (revision 7) |
| II bond operators never removed for J₁ = J₂ = 0 | **fixed** (revision 8) |
| `0/0 = NaN` in the cluster acceptance for C₁ = 1 or C₂ = 1 | **fixed** (revision 9) |
| ZZ formula applied to II bonds in the cluster (found during this round) | **fixed** (revision 9) |
| Non-standard χ_M estimator | **fixed** (revisions 10–12) |
| Importing `data_analysis.py` runs `print_all` (found during this round) | **fixed** (revision 12) |
| A stale `Precomputed/` cache can bring back the old star masks | **fixed** (revision 13) |
| Case 2 untested | **tested**, Sec. 8 |
| `main.jl` / `run_job.jl` path never run | **run**, Sec. 7 (Z and X basis, using HDF5.jl + ArgParse in a throwaway project environment) |
| R2 (missing star-interval move) | **still open:** the move is validated in `../review_scripts/revised*/star_interval_addon.jl` but not yet integrated into `SSE.jl` |
| Slow XX sampling for J₂ ≠ 0 | **open**, Sec. 8.3 |
| L ≥ 3, and the X basis away from J₁ = J₂ | **not tested** |

## 7. Second round: details and verification

Scripts are in `../review_scripts/revised2/`. Run them from that directory;
they read `body.jl` and `apply.jl`, which are extracted from the current
`SSE.jl`.

**Exact cluster acceptance (revision 9).** During the traversal, the cluster
now only records which operators it touches. Once the cluster is complete,
`operator_flip_ratio` returns, for each touched operator, W(after)/W(before).
"Before" means original + previously accepted clusters (`vertices_to_flip`);
"after" also includes this cluster.

- **Pass-through operators (ZZZZ, and ZZ/II if J₂ = 0):** the matrix element
  depends only on the parity of the number of flipped spin lines (legs / 2).
  For ZZZZ, w = C₁ ± B. For ZZ, w = J₁(C₂ ± zz). For II, w = C₄, so the ratio
  is 1.
- **Side-type operators (IIII/XXXX, and ZZ/II/XX if J₂ ≠ 0):** the type is
  toggled iff exactly one side is flipped. A side flip flips both spins of a
  bond, so zz is always the value in the original incoming state.
- A new weight of 0 (possible for C = 1) gives ratio 0 instead of `0/0`.

**Brute-force tests** (`test_cluster_full.jl`): for each cluster, `p_acc`
is compared with W(after)/W(before) computed from the full configuration, and
the flip application is compared with an independent implementation.

| β, J₁, J₂, C₁, C₂ | clusters | ratio mismatches | NaN | zero true ratios | flip application / weights |
|---|---|---|---|---|---|
| 5, 0.1, 0, 2, 2 | 6 048 | 0 | 0 | 0 | 300/300 OK |
| 5, 0.1, 0.1, 2, 2 | 5 043 | 0 | 0 | 0 | 300/300 OK |
| 5, 0.1, 0, **1, 1** | 5 922 | 0 | 0 | 1 017 | 300/300 OK |
| 5, 0.1, 0.1, **1, 1** | 5 474 | 0 | 0 | 103 | 300/300 OK |
| 2, 1, 0, 2, 2 | 977 | 0 | 0 | 0 | 200/200 OK |
| 2, **0, 0**, 2, 2 (pure TC, II bonds) | 2 458 | 0 | 0 | 0 | 300/300 OK, max n = 93 (bounded) |
| 2, **0**, 0.3, 2, 2 | 2 822 | 0 | 0 | 0 | 300/300 OK |

For C = 2 the new acceptance is numerically identical to the old one: a rerun
of the case-1 test with the same seeds reproduces Sec. 4 to every printed
digit.

**Star-interval add-on at J₂ = 0.1** (`test_addon.jl`): 0 ratio mismatches in
17 734 proposals.

**Production path (`main.jl` → `run_job.jl` → `thermalize`/`simulate` →
HDF5 → `data_analysis.py`).** This was run at case 2 in both bases:
thermalize for 5 000 sweeps from `--M 20` (M grew to 211), then simulate for
30 000 sweeps as a warm start with `--use-prev-rng`. The command used:

```
julia --project=<env with HDF5, ArgParse> Julia/main.jl --filename <out>/thermalize_Z.h5 --L 2 --beta 5 \
  --J1 0.1 --J2 0.1 --C1 2 --C2 2 --C3 1 --C4 0 --basis Z --task thermalize --num-sweeps 5000 --state 0 --M 20 --seed 7
julia --project=<env> Julia/main.jl --filename <out>/simulate_Z.h5 ... --task simulate --num-sweeps 30000 \
  --prev-file <out>/thermalize_Z.h5 --use-prev-rng
```

- Writing works (no `DimensionMismatch`). The thermalized configuration read
  back from HDF5 is consistent, periodic and has nonzero weight (86 padding
  slots).
- `measurements/m_n` and `measurements/m_2_n` are present, and
  `data_analysis.compute_all` imports and runs.
- At J₁ = J₂ the X-basis run is identical to the Z-basis run with A ↔ B and
  ZZ ↔ XX exchanged, as it should be.
- Writing the original (wrong) masks into `Precomputed/diagonals_L2.h5` makes
  `run_job.jl` warn and regenerate the correct file.
- Note: the production path does **not** contain the star-interval move yet, so
  its numbers still have the R2 bias (A = 0.933 ± 0.008 in this short run, vs
  0.962).

The ED with the corrected masks gives M_x² = M_z² at J₁ = J₂. This confirms
the electric–magnetic self-duality that the X-basis trick relies on, with the
corrected star geometry: the map h(x,y) → v(x,y−1), v(x,y) → h(x−1,y) sends
stars to plaquettes, plaquettes to stars, and bonds to bonds (checked numerically for L = 2, 3, 4).

## 8. Test of case 2 (L = 2, β = 5, J₁ = J₂ = 0.1)

All SSE runs use 8 independent seeds, 20 000 thermalization + 100 000
measurement sweeps, and fixed M = 300 (driver `case2.jl`). Errors are the
standard error over seeds. χ_M uses the corrected estimator. "Z-loop" is the
diagonal probability of (Z₀Z₁, Z₄Z₆) in the thermal state; these loops are not
conserved for J₂ ≠ 0 (`zloop_ed.py`).

### 8.1 Results

| Observable | ED (corrected model) | Revised package (rev. 1–14) | Revised package + star-interval add-on |
|---|---|---|---|
| E | −8.33792 | −8.277 ± 0.028 (+2.1σ) | **−8.331 ± 0.013** (+0.6σ) |
| C | 0.7638 | 0.77 ± 0.24 | **0.80 ± 0.17** |
| A | 0.96225 | 0.950 ± 0.008 (−1.4σ) | **0.96278 ± 0.00057** (+0.9σ) |
| B | 0.96225 | 0.971 ± 0.006 (+1.4σ) | **0.9613 ± 0.0027** (−0.4σ) |
| ZZ | 0.19997 | 0.215 ± 0.011 | **0.1994 ± 0.0027** (−0.2σ) |
| XX | 0.19997 | 0.156 ± 0.029 | **0.196 ± 0.016** (−0.2σ) |
| M² | 0.29373 | 0.306 ± 0.009 | **0.2938 ± 0.0029** (0.0σ) |
| M⁴ | 0.26351 | 0.278 ± 0.011 | **0.2636 ± 0.0036** (0.0σ) |
| χ_M | 1.61319 | 1.72 ± 0.08 | **1.614 ± 0.023** (0.0σ) |
| P(Z-loop 00) | 0.5537 | 0.576 ± 0.017 | **0.5557 ± 0.0073** |
| P(Z-loop 01 / 10 / 11) | 0.1488 each | 0.141 / 0.142 / 0.142 | **0.148 / 0.147 / 0.149** |
| XX spread over seeds | — | 0.073 … 0.326 | 0.137 … 0.278 |

For comparison, the **original code** (Sec. 6 of `../QMC_REVIEW_FINDINGS.md`)
had XX wandering from 0.05 to 0.87 and M² from 0.008 to 0.52 between 5 000-sweep
blocks.

### 8.2 Reading

1. **Revised package alone:** the averages are within about 2σ, but only
   because the error bars are large. The seed-to-seed spread is big, and in
   10 000-sweep blocks XX sits at about 0.07–0.09 for long stretches and then
   bursts to 0.25–0.55. This is still the student's "not converging" picture
   in a milder form. The sector move helps (the Z-loop weights are no longer
   frozen), but R2 is still present.
2. **With the star-interval add-on:** every observable agrees with ED within
   1σ, and A and M² have small error bars (0.06% and 1%). So for case 2 the
   remaining R2 bias is also cured by this move.
3. **XX is still the slowest channel** (relative error 8%, blocks between 0.09
   and 0.6). The error bars are honest (the seed spread matches the in-run
   block spread), so the result is correct but expensive.

### 8.3 Most probable reason for the slow XX channel, and proposed next step

This is the bond analogue of R2. For J₂ ≠ 0, bond operators (ZZ/XX) are
side-type in the cluster: a line reaching them always takes both legs of that
side. A single spin line therefore cannot pass through a bond operator, and a
pair of XX operators on the same bond can only be created or annihilated by a
cluster that happens to connect them directly. The proposed **bond-interval
move** follows the star-interval construction (not implemented yet):

- pick a bond b and an ordered pair of bond operators (ZZ or XX) on b;
- toggle both (ZZ ↔ XX);
- flip the two spins of b on the time slices between them.

The weight ratio comes from the toggles (J₂ vs J₁(C₂ + zz)). Between the two
operators, it also includes the ZZZZ operators of the two plaquettes that
contain exactly one spin of b, and the other ZZ bonds with odd overlap. In
other words it creates a pair of fluxes for the interval, as a physical XX
pair does.

## 9. Recommended next steps (after round 2; see Sec. 10 for what was done)

1. ~~Integrate the star-interval move into `SSE.jl`.~~ Done (revision 15).
2. ~~Implement and test the bond-interval move.~~ Done (revision 16).
3. Test the X basis away from J₁ = J₂, and L = 3 (third case of `info.txt`). L = 3 and the case-2 tension: see Sec. 13.
4. ~~Make the interval moves O(M) per sweep (Sec. 10.5).~~ Done (Sec. 11).

## 10. Third round: interval moves integrated (2026-10-24)

| # | File | Change |
|---|---|---|
| 15 | `Julia/SSE.jl`: new `interval_ratio`, `interval_apply!`, `interval_update!`, `star_interval_update!`, called at the end of `qmc_sweep` | The star-interval move (validated earlier as an add-on) is now part of every sweep: one proposal per star. |
| 16 | `Julia/SSE.jl`: new `bond_interval_update!`, called at the end of `qmc_sweep` | **New bond-interval move** (for J₂ ≠ 0): one proposal per NN bond. Pick an ordered pair of bond operators (ZZ/II or XX) on the bond, toggle both (ZZ/II ↔ XX), and flip the bond's 2 spins on the slices between them. Inside the interval, the two plaquettes containing exactly one spin of the bond flip their B_p (a flux pair, like a physical XX pair). Diagonal operators with odd overlap enter the ratio; ZZ endpoints are evaluated on their original incoming state. |

Both moves share one kernel. The comment block above `interval_ratio` in
`SSE.jl` gives the consistency, weight-ratio and proposal-symmetry argument.
The star-interval move is unchanged from the add-on: case 1 run through the
package alone reproduces the add-on results of Sec. 4 to every printed digit.

### 10.1 Brute-force test of the integrated moves (`../review_scripts/revised3/test_interval_moves.jl`)

For random (a, c) proposals, `interval_ratio` is compared with W(after)/W(before)
computed from the full configuration. The configuration is also checked after
every full `qmc_sweep`.

| β, J₁, J₂, C₁, C₂ | star proposals (mismatches) | bond proposals (mismatches) | invalid configs |
|---|---|---|---|
| 5, 0.1, 0.1, 2, 2 | 8 498 (0) | 7 212 (0) | 0/800 |
| 5, 0.1, 0.1, 1, 1 | 8 498 (0) | 4 235 (0), 4 041 with ratio 0 | 0/800 |
| 5, 0.1, 0, 2, 2 | 8 411 (0) | — (off for J₂ = 0) | 0/800 |
| 2, 0, 0.3, 2, 2 (II ↔ XX) | 6 283 (0) | 17 731 (0) | 0/800 |
| 1, 1, 1, 2, 2 | 2 575 (0) | 21 287 (0) | 0/800 |
| 2, 0, 0, 2, 2 | 6 232 (0) | — | 0/800 |

### 10.2 Exact test of the complete transition kernel (truncated ensemble)

This test is independent of how fast the chain mixes. With a fixed small cutoff
M, the SSE samples exactly
Z_M = Σ_{n≤M} βⁿ/n! Tr Tⁿ, with T = C − H, the sum of all operator weights.
Its expectations follow exactly from matrix powers (`truncated_exact.py`). A
short chain mixes fast, so all moves together (diagonal, cluster, column flips,
loop move, star- and bond-interval) can be tested very precisely. The driver
is `truncated_sse.jl` (12 seeds × 3·10⁶ sweeps, after 2·10⁴ thermalization).

| quantity | β = 1, J₁ = J₂ = 0.5, M = 14: SSE / exact (pull) | β = 2, J₁ = J₂ = 0.3, M = 30: SSE / exact (pull) |
|---|---|---|
| ⟨n⟩ | 13.42869(16) / 13.42869 (0.0σ) | 29.17530(35) / 29.17522 (+0.2σ) |
| ⟨n_XX⟩ | 0.8391(22) / 0.8360 (+1.5σ) | 1.8264(77) / 1.8373 (−1.4σ) |
| ⟨n_XXXX⟩ | 0.6300(15) / 0.6302 (−0.1σ) | 2.7476(34) / 2.7497 (−0.6σ) |
| P(Z-loop 00) | 0.41763(25) / 0.41759 (+0.2σ) | 0.57801(67) / 0.57727 (+1.1σ) |
| P(Z-loop 01) | 0.20014(8) / 0.20005 (+1.2σ) | 0.14198(22) / 0.14208 (−0.5σ) |
| P(Z-loop 10) | 0.19987(10) / 0.20005 (−1.9σ) | 0.14186(26) / 0.14208 (−0.9σ) |
| P(Z-loop 11) | 0.18237(14) / 0.18231 (+0.4σ) | 0.13815(23) / 0.13857 (−1.8σ) |
| ⟨m²⟩ (slot 1) | 0.36150(26) / 0.36160 (−0.4σ) | 0.46899(61) / 0.46833 (+1.1σ) |
| P(n = M) | 0.62516(10) / 0.62513 (+0.3σ) | 0.53711(13) / 0.53701 (+0.8σ) |

A third point (β = 1.5, J₁ = J₂ = 0.1, M = 14, 8 seeds × 10⁶) also agrees
within 1.3σ everywhere. **Conclusion: the full kernel samples the SSE weight
exactly** (to about 10⁻⁴ in the sector probabilities).

### 10.3 Case 1 and case 2 through the package alone (`qmc_sweep` only, no add-on)

- **Case 1** (8 seeds × 10⁵ sweeps): identical to Sec. 4, last column. Every
  observable is within 1.5σ of ED.
- **Case 2** (28 seeds × 10⁵ sweeps, `case2_pkg.jl`):

| Observable | ED | Package (star + bond interval) | Previous: star interval only (8 seeds) |
|---|---|---|---|
| E | −8.33792 | −8.3357 ± 0.0039 (+0.6σ) | −8.331 ± 0.013 |
| A | 0.96225 | 0.96231 ± 0.00042 (+0.1σ) | 0.96278 ± 0.00057 |
| B | 0.96225 | 0.96331 ± 0.00093 (+1.1σ) | 0.9613 ± 0.0027 |
| ZZ | 0.19997 | 0.20159 ± 0.00079 (+2.0σ) | 0.1994 ± 0.0027 |
| XX | 0.19997 | 0.1939 ± 0.0042 (−1.4σ) | 0.196 ± 0.016 |
| M² | 0.29373 | 0.29519 ± 0.00073 (+2.0σ) | 0.2938 ± 0.0029 |
| M⁴ | 0.26351 | 0.26525 ± 0.00087 (+2.0σ) | 0.2636 ± 0.0036 |
| χ_M | 1.61319 | 1.6222 ± 0.0060 (+1.5σ) | 1.614 ± 0.023 |
| P(Z-loop 00 / 01 / 10 / 11) | 0.5537 / 0.1488 ×3 | 0.5571(18) / 0.1466(6) / 0.1486(8) / 0.1477(7) | — |

Per 8 seeds, the XX error dropped from 0.016 to about 0.005 (XX seed spread
0.17–0.21 instead of 0.14–0.28), and the other error bars dropped by a factor
2–3. **Most of the improvement is from the bond-interval move.**

**Remaining tension:** with 28 seeds, the sector-00 weight and the observables
correlated with it (ZZ, M², M⁴) sit about 2σ high, and P(Z-loop 01) is 3.7σ
low. Sec. 10.2 rules out a detailed-balance error, because the same moves
sample the exact truncated ensemble. The most likely explanation is that the
slowest mode (XX / Z-loop sector changes at β = 5) has an autocorrelation time
that is not small compared with 10⁵ sweeps. The per-seed averages then are not
fully independent of the starting configuration (each run starts from a single
random basis state), and the error bars are somewhat optimistic. To settle
it: longer runs (10⁶ sweeps per seed) and an integrated autocorrelation time
for the Z-loop sector indicator.

### 10.4 Cost of the interval moves (`bench_interval.jl`)

| L | β | full sweep | interval moves alone | share |
|---|---|---|---|---|
| 2 | 2 / 5 / 10 | 0.063 / 0.161 / 0.327 ms | 0.002 / 0.006 / 0.016 ms | 3.4 / 3.6 / 4.9% |
| 3 | 2 / 5 | 0.139 / 0.342 ms | 0.009 / 0.026 ms | 6.5 / 7.6% |
| 4 | 2 / 5 | 0.254 / 0.709 ms | 0.024 / 0.087 ms | 9.5 / 12.3% |

The share grows roughly ∝ L², as expected for O(N_masks · M) = O(βL⁴) against
O(M) = O(βL²). Negligible at L = 2, but dominant for large L.

### 10.5 Planned optimization (implemented in Sec. 11)

This follows a review comment; the comment's analysis is correct.

1. Build, once per interval stage, a list of operator positions for every mask
   (star or bond) and for every spin, in one O(M) pass. Positions and masks
   do not change during the interval moves (toggling keeps both), so the lists
   stay valid. Finding a pair is then O(1).
2. Compute the ratio from the operators on the mask's spins only (2 or 4
   per-spin lists, each operator counted once), not from every slice in the
   interval.
3. **Point the comment does not mention:** `interval_apply!` also flips the
   per-slice `states` array over the whole interval, which is O(interval
   length) per accepted move. To reach O(βL²), keep the current eigenvalue
   (B_p, or Z_iZ_j of the bond) of each operator in a per-operator array. An
   accepted move then updates only the operators it walked, and `states` is
   rebuilt once per sweep by propagating from `states[1]`, which is flipped by
   S if an accepted interval wraps through slot 1.
4. Cache the lattice lists (stars, bonds, plaquettes, dual loops) per L instead
   of rebuilding them every sweep.
5. Optional: choose c as the next operator on the same mask after a
   (cyclically). This is still symmetric, because toggling does not move
   operators. It gives short intervals and higher acceptance, but changes the
   Markov chain, so it must be judged on accuracy per CPU time.

Items 1–4 keep the proposal distribution and the order of random numbers
unchanged, so the results should reproduce the current ones (up to rounding
in the ratio products). That gives a direct regression test.

## 11. Fourth round: optimized interval stage and the nearby-pair option (2026-10-24)

| # | File | Change |
|---|---|---|
| 17 | `Julia/SSE.jl`: new `interval_stage!` (replaces `interval_ratio`, `interval_apply!`, `interval_update!`, `star_interval_update!`, `bond_interval_update!`) | **Optimized interval stage.** (i) One O(M) pass builds the operator list of every mask (star or bond), the sorted list of operators on every spin, and each operator's current eigenvalue `sgn` (B_p for ZZZZ, Z_iZ_j for ZZ/XX). (ii) A proposal walks only the operators on the spins of S inside the interval, found by binary search in the per-spin lists. Each operator is counted once (at its lowest spin in S). (iii) An accepted move flips only the `sgn` of the operators it walked, plus the endpoint types. (iv) `states` is rebuilt once at the end of the stage by propagating from `states[1]`, which is first flipped by S for every accepted interval that wraps through slot 1. |
| 18 | `Julia/SSE.jl`: `interval_stage!(...; pairs)`, `qmc_sweep(...; interval_pairs)` | **Pair choice option.** `:uniform` (**default**): a, c independently uniform among the operators on S, i.e. the original proposal. `:nearby`: a uniform, c = the next operator on S after a (cyclically). Both are symmetric, because toggling moves no operator. |
| 19 | `Julia/SSE.jl`: new `lattice_lists(L)` with `LATTICE_CACHE` | Stars, bonds and dual loops are cached per L instead of being rebuilt every sweep (used by the loop move and the interval stage). |

The production path (`thermalize`, `simulate`) calls `qmc_sweep` without the
keyword, so it uses the default `:uniform`.

### 11.1 Minimal tests (`../review_scripts/revised4/`)

1. **Regression (`regress.jl`):** the optimized code with `:uniform`
   reproduces the pre-optimization chain (`SSE_before_optimization.jl`)
   **bit for bit**: identical operator counts and an identical hash of the
   final configuration at (β, J₁, J₂, M) = (5, 0.1, 0.1, 300),
   (5, 0.1, 0, 300) and (1, 0.5, 0.5, 14). The default `qmc_sweep(...)`
   equals `interval_pairs = :uniform`.
2. **Exactness of `:nearby` against the truncated ensemble (`truncated_sse.jl`,
   `truncated_exact.py`):**
   - β = 1, J₁ = J₂ = 0.5, M = 14, 12 seeds × 3·10⁶ sweeps: all 9 quantities
     within 1.7σ, e.g. P(Z-loop 00) = 0.41756(10) vs 0.41759, ⟨n⟩ =
     13.42899(17) vs 13.42869.
   - β = 2, J₁ = J₂ = 0.3, M = 30, 6 seeds × 10⁶: all within 1.7σ.
   - (A first 6-seed run at the first point showed a 3.9σ pull. It
     disappeared with more statistics, as for `:uniform` in Sec. 10.2.)
3. **Cost (`bench_opt.jl`, J₁ = J₂ = 0.1), interval stage per sweep:**

| L, β | old | optimized `:uniform` | optimized `:nearby` | stage share of an optimized sweep |
|---|---|---|---|---|
| 2, 5 | 0.009 ms | 0.018 ms | 0.015 ms | 8% |
| 2, 10 | 0.021 ms | 0.050 ms | 0.030 ms | 8% |
| 3, 5 | 0.041 ms | 0.048 ms | 0.032 ms | 8% |
| 4, 5 | 0.178 ms | 0.097 ms | 0.075 ms | 10% |
| 5, 5 | 0.473 ms | 0.158 ms | 0.121 ms | 9% |

   The optimized stage is now a constant ~8–10% of a sweep (O(M) scaling),
   while the old one grew ∝ L² relative to the sweep. At L = 2 the optimized
   version is slightly slower in absolute terms (fixed setup cost of the lists
   and of the rebuild), which is irrelevant.
4. **Efficiency of `:nearby` vs `:uniform` in case 2** (same 8 seeds × 10⁵
   sweeps): XX = 0.1843 ± 0.0089 (`:nearby`) vs 0.1939 ± 0.0048 (`:uniform`);
   P(Z-loop 00) error 0.0040 vs 0.0028; A and E comparable. Short intervals
   decorrelate the slow XX / sector mode less well than uniformly chosen
   pairs. Since the cost difference is at most about 30%, **`:uniform` gives
   better accuracy per CPU time here and is the default.** `:nearby` remains
   available (it may pay off at large L or β, where long intervals are rarely
   accepted).

### 11.2 Notes

- The spin configuration is stored as one `Int64` bitmask, so the code (the
  original as well) works only for N = 2L² ≤ 63, i.e. **L ≤ 5**. At L = 6 the
  benchmark fails inside the cluster construction. Larger L needs a
  multi-word spin representation.
- Still open: the case-2 sector tension of Sec. 10.3 (longer runs and an
  autocorrelation time), the X basis away from J₁ = J₂, and L = 3.

## 12. Housekeeping after round 4 (2026-10-24)

- **Line endings of `Julia/SSE.jl`.** Every file of the original package uses
  Windows (CRLF) line endings. One scripted edit in round 4 had silently
  converted `SSE.jl` to Unix (LF), so `git diff` against the original showed
  every line as changed. The file was converted back to CRLF. Check: it loads
  and reproduces the same chain (same final-configuration hash as before the
  conversion). No code changed.
- **New `REVISION_SUMMARY.md`.** A reader-oriented summary of all changes
  relative to the original package, with the reason for each change. This log
  stays the detailed, chronological record. The summary is the place to start.
- No further code changes. Open after round 4: L = 3, the case-2 sector
  tension (Sec. 10.3), and the X basis away from J₁ = J₂. Sec. 13 covers the
  first two.

## 13. Fifth round: tests of L = 3 and of the case-2 tension (2026-10-24, no code changes)

Scripts and raw output are in `../review_scripts/revised5/`:

- `ed_sym.jl`: an independent exact diagonalization. It takes the lattice masks
  from `SSE.jl`. Everything else is new: the symmetry basis, matrix elements,
  observables, and the Kubo χ_M. H is block-diagonalized with translations
  (L × L momenta), the global spin flip ∏X and the Z parity ∏Z. Each block is
  then diagonalized exactly. At L = 3 there are 36 blocks of 7280–7296 states,
  whose dimensions sum to exactly 2¹⁸. The k and −k blocks are complex
  conjugates, so 20 diagonalizations are enough. Run time: 40 min. Diagonal
  observables that are not symmetric, such as the Z-loop sector projector at
  L = 2, are averaged over the symmetry orbit, which is exact because ρ commutes
  with the group. **Validation:** at L = 2 it reproduces the earlier ED of case
  1 and case 2 to every printed digit. At L = 3 and J₁ = J₂ it gives exactly
  A = B and ZZ = XX, as self-duality requires.
- `mc_run.jl`: a package-only run (`qmc_sweep` from this `SSE.jl`) for any
  L ≤ 5. The constants are those of `make_jobs.py`. Thermalization grows the
  cutoff like `thermalize.jl` (M ≥ 4n/3), and the measurement phase keeps M
  fixed. It records exact running sums and 10-sweep bin averages.
- `analyze.py`: means over seeds with the standard error from the spread between
  seeds, pulls against ED, integrated autocorrelation times τ_int from the bin
  series (Sokal windowing), and a drift check (first 10% of the measurement
  phase vs. the rest).
- `case2_pkg_old_protocol.jl`: the Sec. 10.3 driver, unchanged (fixed M = 300,
  2·10⁴ thermalization sweeps). It records all four Z-loop sector probabilities.
- `typicality.jl`: a thermal-typicality estimate on the full sparse H, checked
  at L = 2. It was not needed in the end, because the exact L = 3 result was
  already available.

### 13.1 L = 3, β = 1, J₁ = J₂ = 1 (third case of `info.txt`)

16 seeds × (10⁵ thermalization + 10⁶ measurement) sweeps, `mc_run.jl`:

| Observable | ED | SSE | pull |
|---|---|---|---|
| E | −49.010475 | −49.0064 ± 0.0093 | +0.4σ |
| C | 0.493904 | 0.446 ± 0.055 | −0.9σ |
| A | 0.545833 | 0.5485 ± 0.0041 | +0.7σ |
| B | 0.545833 | 0.5425 ± 0.0040 | −0.8σ |
| ZZ | 0.544243 | 0.5410 ± 0.0040 | −0.8σ |
| XX | 0.544243 | 0.5476 ± 0.0042 | +0.8σ |
| M² | 0.504314 | 0.5010 ± 0.0041 | −0.8σ |
| M⁴ | 0.421496 | 0.4180 ± 0.0043 | −0.8σ |
| χ_M | 7.890887 | 7.82 ± 0.08 | −0.8σ |

**The revised code is correct at L = 3.** It mixes slowly, though:

| | n | n_A, n_B | n_ZZ, n_XX, M², M⁴, χ_M |
|---|---|---|---|
| τ_int (sweeps, median over seeds) | ≈ 12 | ≈ 620–630 | ≈ 800–850 |

The slow mode is an anticorrelated exchange of weight between the X-type
terms (A, XX) and the Z-type terms (B, ZZ, M²). Seed means range from
A = 0.586, B = 0.506 (seed 309, "X-heavy") to A = 0.513, B = 0.577 (seed 303,
"Z-heavy"), while A = B exactly. At J₁ = J₂ = 1 both Ising perturbations are
as strong as the stabilizer terms, so Z-type and X-type ordering compete, and
the SSE moves between the two only slowly. There is no sign of a bias: no
drift after thermalization (all |pull| ≤ 0.5σ), and the spread between seeds
is about 1.2× what τ_int predicts.

**Reading of the original report** ("does not appear to be converging, the M²
plot appears weird"): the original code at L = 3 simulated a different model,
whose "stars" anticommute with the plaquettes (Sec. 1, revision 1). Even with
the correct model, though, a run at this point has τ_int ≈ 10³ sweeps, so an
M² time series shows long plateaus, and runs of 10⁴–10⁵ sweeps look
unconverged. **Advice:** at this point use ≥ 10⁶ sweeps per seed and several
seeds, and monitor n_A − n_B or M² in blocks of ≥ 10⁴ sweeps.

### 13.2 Case 2 (L = 2, β = 5, J₁ = J₂ = 0.1): the Sec. 10.3 tension

Two independent sets of 16 seeds × 10⁶ measurement sweeps each, 10× longer per
seed than in Sec. 10.3. Run 1 used `mc_run.jl` (grown cutoff, 10⁵
thermalization sweeps). Run 2 used `case2_pkg_old_protocol.jl` (the Sec. 10.3
protocol). The exact P(Z-loop) values come from `ed_sym.jl`.

| Observable | ED | Run 1 | Run 2 | Combined (32 seeds) | pull |
|---|---|---|---|---|---|
| E | −8.337921 | −8.33829(137) | −8.33751(152) | −8.33794(102) | −0.0σ |
| A | 0.962251 | 0.962146(123) | 0.962187(182) | 0.962159(102) | −0.9σ |
| B | 0.962251 | 0.961737(359) | 0.963063(320) | 0.962476(239) | +0.9σ |
| ZZ | 0.199973 | 0.199959(328) | 0.200092(434) | 0.200008(262) | +0.1σ |
| XX | 0.199973 | 0.20183(168) | 0.19778(160) | 0.19971(116) | −0.2σ |
| M² | 0.293732 | 0.293722(295) | 0.293906(306) | 0.293810(212) | +0.4σ |
| M⁴ | 0.263513 | 0.263417(357) | 0.263766(363) | 0.263588(254) | +0.3σ |
| χ_M | 1.613185 | 1.61415(245) | 1.61475(260) | 1.61443(178) | +0.7σ |
| P(Z-loop 00) | 0.553697 | 0.553594(742) | 0.554209(719) | 0.553911(517) | +0.4σ |
| P(Z-loop 01) | 0.148768 | — | 0.148703(308) | | −0.2σ |
| P(Z-loop 10) | 0.148768 | — | 0.148273(301) | | −1.6σ |
| P(Z-loop 11) | 0.148768 | — | 0.148815(259) | | +0.2σ |

**The Sec. 10.3 tension was a statistical fluctuation.** Evidence:

- **Agreement with ED:** with 10× the statistics, the quantities that were 2σ
  high (P(Z-loop 00), ZZ, M², M⁴) agree with ED within 0.7σ, and P(Z-loop 01),
  which was −3.7σ, is now −0.2σ.
- **Error bars:** the error of P(Z-loop 00) dropped from 0.0018 (28 × 10⁵
  sweeps) to 0.00074 (16 × 10⁶), as √(sweeps) predicts, so the earlier error
  bars were honest.
- **The suspected mechanism is ruled out:** the Z-loop sector is not slow. Its
  τ_int is ≈ 15 sweeps. The slowest mode at this point is n_XX with τ_int ≈ 460
  sweeps (max 630). Every other count has τ_int ≈ 8–16 sweeps.
- **No drift** after thermalization.

The B difference between the two runs (2.8σ) goes with an opposite XX
difference (+1.8σ). An XX operator anticommutes with the plaquettes it
touches, so B shares a small part of the slow XX mode, and with 16 seeds its
error bar is itself uncertain by about ±18%. Combined, B is +0.9σ from ED.

### 13.3 Status

- L = 3: done (13.1). The code is correct, and τ_int ≈ 10³ sweeps at β = 1,
  J₁ = J₂ = 1.
- Case-2 tension: resolved (13.2).
- Still untested: the X basis away from J₁ = J₂.
