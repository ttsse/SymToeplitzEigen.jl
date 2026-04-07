# SymToeplitzEigen

<div align="center">

| **Build Status**                                                                                |
|:-----------------------------------------------------------------------------------------------:|
| [![CI](https://github.com/ttsse/SymToeplitzEigen.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ttsse/SymToeplitzEigen.jl/actions/workflows/CI.yml) |

</div>





## Installation

```julia
pkg> add https://github.com/ttsse/SymToeplitzEigen.jl
```
or
```julia
julia> using Pkg
julia> Pkg.add("https://github.com/ttsse/SymToeplitzEigen.jl")
```
___
## Usage

This package implements solving eigenvalues and eigenvectors of symmetric Toeplitz-like matrices in low precision and then using Iterative Refinement (inspired by [IterativeRefinement.jl](https://github.com/RalphAS/IterativeRefinement.jl)) in parallel to improve the accuracy.

### Exported functions:

```julia
vals, vecs = EigenRef(A)
```

where `A` is a symmetric Toeplitz-like matrix.

```julia
vals, vecs = EigenRef(n, vc)
```

where `n` is the size of the matrix, `vc` is the first column of symmetric Toeplitz matrix.

Also supported:

```julia
using SparseArrays
A = spdiagm(-1 => fill(-1.0, 99), 0 => fill(2.0, 100), 1 => fill(-1.0, 99))
vals, vecs = EigenRef(A)
```

If optional weak dependencies are loaded, `EigenRef` also dispatches directly on:

- `ToeplitzMatrices.SymmetricToeplitz` / symmetric `ToeplitzMatrices.Toeplitz`
- `BandedMatrices.AbstractBandedMatrix`

<!-- #Refinement_precision :: Integer = 256, Max_iter :: Integer = 10, tol_fact :: Integer = 1 -->

### Optional keyword arguemnts:

* `Low_pres_type :: Type = Float32`

    The data type to do the initial low precision computation. Should either be `Float32` or `Float64`. Not used if given matrix is already one of these types.
  

<br />

* `Refinement_precision :: Integer = 256`
    
    The precision up to which the eigenvalues and eigenvectors should be computed. Automatically computed if the given symmetric Toeplitz-like matrix is of type `BigFloat`.

<br />

* `Max_iter :: Integer = 100`

    The maximum number of iterations of Iterative Refinement to perform.

    <br />

* `tol_fact :: Integer = 1`

    The tolerance factor up to which the error is computed, that is the error estimate is computed up to `tol_fact*eps(BigFloat)` where `BigFloat` has precision `Refinement_precision`

    <br />

* `solve_mode :: Symbol = :fast`

    Solve policy used for the correction system inside refinement.
    - `:fast` uses `Float64` factorization/solves and promotes corrections to `BigFloat`
    - `:robust` performs factorization/solves directly in `BigFloat`
    - `:adaptive` starts in `:fast` and promotes a pair to `:robust` when progress stalls

    <br />

* `stall_ratio :: Float64 = 0.9`

    For `solve_mode = :adaptive`, a pair is considered stalled when the relative residual reduction factor is greater than or equal to `stall_ratio` for consecutive iterations.

    <br />

* `stall_iters :: Integer = 3`

    Number of consecutive stalled iterations before an adaptive pair is promoted from `:fast` to `:robust` mode.

    <br />

* `return_status :: Bool = false`

    If `true`, `EigenRef` returns a third value `status` with per-eigenpair metadata: `mode_used`, `iterations_used`, `converged`, `promoted_iteration`, `precision_bits_used`, and `precision_escalated`.

    <br />

* `show_progress :: Bool = false`

    If `true`, a live progress indicator is displayed while eigenpairs are refined.
    If `ProgressMeter.jl` is already loaded in the active Julia session, it is
    used automatically; otherwise a built-in text progress bar is shown.

    <br />

* `adaptive_precision_escalation :: Bool = false`

    If `true` and `solve_mode = :adaptive`, unconverged pairs after base-precision refinement are retried in robust mode at higher precision. This keeps extra memory and compute focused on hard pairs only.

    <br />

* `escalation_precision :: Integer = 0`

    Target precision in bits for escalated pairs. If set to `0`, the package uses `2*Refinement_precision`.

    <br />

* `escalation_extra_iter :: Integer = 10`

    Maximum number of additional robust iterations used for escalated pairs.

    <br />

* `toeplitz_kernel :: Symbol = :auto`

    Selects how residual/factorization base matrices are formed in refinement.
    - `:auto` selects `:structured` only when matrix size is above `toeplitz_auto_threshold` and the matrix is Toeplitz (or `EigenRef(n, vc)` provides `vc`)
    - `:dense` uses the dense matrix path
    - `:structured` forces symmetric Toeplitz kernels (requires Toeplitz input or provided `vc`)

    <br />

* `toeplitz_auto_threshold :: Integer = 96`

    Minimum matrix size where `toeplitz_kernel = :auto` considers switching to the structured Toeplitz path.

    <br />

* `reuse_toeplitz_cache :: Bool = true` (only for `EigenRef(n, vc)`)

    Reuses a cached Toeplitz matrix when repeated calls use the same `n` and `vc`, avoiding matrix reconstruction overhead.

    <br />

* `sparse_solver :: Symbol = :auto` (sparse inputs only)

    Selects refinement correction backend for sparse matrices.
    - `:auto` uses sparse solves when matrix density is low and `solve_mode = :fast`
    - `:sparse` forces sparse correction solves
    - `:dense` forces dense fallback refinement

    <br />

### Nonsymmetric API

```julia
vals, right_vecs, left_vecs = EigenRefNonSym(A)
```

or with explicit initial guesses:

```julia
vals, right_vecs, left_vecs = EigenRefNonSym(A, vals0, vecs0)
```

Use `return_status = true` to also return diagnostic status metadata per pair,
including mode used, convergence, right/left residuals, backward/biorthogonality
errors, condition proxies, clustering, and warnings for near-defective or
unstable cases.

`show_progress = true` is also supported for `EigenRefNonSym`.

### Weak-Dependency Matrix Extensions

The package uses Julia extensions (weak dependencies) for optional matrix ecosystems.
To activate these methods, load the corresponding package in your session:

```julia
using ToeplitzMatrices
using BandedMatrices
```

Then calls like `EigenRef(SymmetricToeplitz(...))`, `EigenRef(BandedMatrix(...))`,
are enabled.

### Safety Gates and Known Failure Regimes

* Symmetric path (`EigenRef`):

    - `toeplitz_kernel = :structured` requires a symmetric Toeplitz input (or explicit Toeplitz first-column data). If this is not satisfied, the structured path is rejected.
    - `toeplitz_kernel = :auto` safely falls back to dense mode when Toeplitz structure is not detected.
    - `solve_mode = :adaptive` can promote stalled pairs to robust solves and, if enabled, to higher precision; this improves reliability but may increase runtime.

* Nonsymmetric path (`EigenRefNonSym`):

    - Guaranteed convergence is limited to diagonalizable, non-clustered target pairs.
    - Clustered or near-defective spectra trigger Schur-block fallback and warning metadata.
    - Strongly non-normal matrices are best-effort; monitor `condition_proxy`, residuals, and warnings.

<br />
___
## Example

```julia
julia> using LinearAlgebra

julia> using SymToeplitzEigen

julia> v = rand(10);

julia> A = SymToeplitzEigen.toeplitz(100, v, v);

julia> nvals, nvecs = EigenRef(A);

julia> maximum([norm((A - nvals[kk]I)*nvecs[:,kk]) for kk in 1:100])
4.631891758490566668625685607982804238714754608769043250856766307722715870331231e-75
```

___
## Testing
After installation, run
```julia
pkg> test SymToeplitzEigen
```
or
```julia
julia> using Pkg
julia> Pkg.test("SymToeplitzEigen")
```

## Performance Baseline

To track runtime, allocations, and residual quality during optimization work, run:

```julia
julia --project -e 'include("test/perf_baseline.jl"); run_suite()'
```

For quicker local checks, pass smaller cases:

```julia
julia --project -e 'include("test/perf_baseline.jl"); run_suite(ns=(1000,), precs=(256,), repeats=1)'
```

To benchmark only the dense path, set kernels explicitly:

```julia
julia --project -e 'include("test/perf_baseline.jl"); run_suite(ns=(1000,), precs=(256,), kernels=(:dense,), repeats=1)'
```

## Plan Benchmark Suite (Steps 1-34)

To benchmark the major optimization tracks from the implementation plan
(allocation/threading, solve policies, Toeplitz kernels/cache,
nonsymmetric path, and thread scalability), run:

```julia
julia --project test/plan_benchmarks.jl --preset=quick --threads=1,2
```

For a larger run closer to release-scale checks:

```julia
julia --project test/plan_benchmarks.jl --preset=full --threads=1,2,4 --repeats=2
```

The suite writes a CSV report to `test/benchmark_results/` by default.
You can set an explicit output path:

```julia
julia --project test/plan_benchmarks.jl --preset=quick --threads=1,2 --out=test/benchmark_results/my_run.csv
```

## Phase 7 Release Gates

Run the automated Phase 7 verification gates (steps 7.1 to 7.3):

```julia
julia --project test/phase7_release_gates.jl
```

This checks:

- runtime and peak-RSS improvement at large `n` using baseline vs optimized symmetric configurations,
- residual/backward quality plus reproducibility across thread counts,
- nonsymmetric right/left residual, biorthogonality, and condition/warning reporting requirements.

For a quicker local smoke run:

```julia
julia --project test/phase7_release_gates.jl --n-perf=1000 --n-quality=400 --max-iter=20 --threads=1,2
```
