# SymToeplitzEigen

<div align="center">

| **Build Status**                                                                                |
|:-----------------------------------------------------------------------------------------------:|
| [![CI](https://github.com/ttsse/SymToeplitzEigen.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ttsse/SymToeplitzEigen.jl/actions/workflows/CI.yml) |

</div>





## Installation

In Julia 1.9
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
