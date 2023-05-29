# ToeplitzCUDA_Eigen_Refinement

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ttsse.github.io/ToeplitzCUDA_Eigen_Refinement.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ttsse.github.io/ToeplitzCUDA_Eigen_Refinement.jl/dev/)
[![Build Status](https://github.com/ttsse/ToeplitzCUDA_Eigen_Refinement.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ttsse/ToeplitzCUDA_Eigen_Refinement.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ttsse/ToeplitzCUDA_Eigen_Refinement.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ttsse/ToeplitzCUDA_Eigen_Refinement.jl) -->

This package implements solving eigenvalues of symmetric Toeplitz matrices in low precision on the GPU (see [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl)) and then using Iterative Refinement (see [IterativeRefinement.jl](https://github.com/RalphAS/IterativeRefinement.jl)) in parallel to imrpove the accuracy.

Exports the function:

```julia
vals, vecs = compute_Toeplitz_Eigen(n, vc)
```

where `n` is the size of the matrix, `vc` is the first column.
