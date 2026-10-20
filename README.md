# Toric Code QMC

Quantum Monte Carlo (QMC) Stochastic Series Expansion (SSE) simulations.

## Structure

- `Julia/` SSE simulation code
  - `main.jl` Accepts command line inputs for a job, passes to run_job
  - `run_job.jl` Given the output filename and job description, retrieve / construct the necessary inputs and run thermalize/simulate.jl 
  - `thermalize.jl` Runs SSE.jl for a given number of sweeps, and pass the raw data and metadata to write_file.jl
  - `simulate.jl` Same as above, except the operator count cutoff is fixed
  - `SSE.jl` Quantum Monte Carlo
  - `write_file.jl` Records the raw data and meta data in a h5 file
- `Python/` job generation and data analysis
  - `data_analysis` Given the filename, compute and/or print all the observables with error bars
  - `thermalization_analysis` Create thermalization figures for all thermalization runs in RawData (unless figures already present)
  - `jobs` Hardcoded list of jobs
  - `make_jobs` Given a basic job type and parameters, construct the full command required by Julia
  - `run_local` Run all the jobs locally
  - `run_pbs` Construct a PBS file and call qsub
<!-- These are not relevant for the final project, but were used during the validation stage 
  - `exact_diagonalization` L=2 ED
  - `sse_vs_ed` Compute the results of SSE and ED, store in a docx file
-->
- `RawData/` raw simulation output
- `Analysis/` figures, summaries 