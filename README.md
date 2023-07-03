# ToeplitzCUDA_Eigen_Refinement

This package implements solving eigenvalues of symmetric Toeplitz matrices in low precision on the GPU (see [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl)) and then using Iterative Refinement (see [IterativeRefinement.jl](https://github.com/RalphAS/IterativeRefinement.jl)) in parallel to improve the accuracy.

Exports the function:

```julia
vals, vecs = compute_Toeplitz_Eigen(n, vc)
```

where `n` is the size of the matrix, `vc` is the first column.
