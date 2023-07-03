module ToeplitzCUDA_Eigen_Refinement


using LinearAlgebra

include("refinement.jl")
include("Cudacompute.jl")
include("helper.jl")

export compute_Toeplitz_Eigen

"""
    compute_Toeplitz_Eigen(n, vc; CUDA_type = Float32, Refinement_precision = 256, Refinement_Iterations = :single)

Computes the eigenvalues (and vectors) of the symmetric Toeplitz matrix of size n x n with the first column given by vector vc. The eigenvalues (and vectors) are first computed in "low" precision (specified using CUDA_type with values of either Float32 or Float64) on the GPU. The calculated values are then improved using an Iterative Refinement procedure in parallel.
"""
function compute_Toeplitz_Eigen(n :: Integer, vc :: Array; CUDA_type :: Type = Float32, Refinement_precision :: Integer = 256, Max_iter :: Integer = 10, tol_fact :: Integer = 1)

    Tn = toeplitz(n, CUDA_type.(vc), CUDA_type.(vc))

    CUDAvals, CUDAvecs = Compute_CUDA_Eigen(Tn)

    return Refinement(Tn, CUDAvals, CUDAvecs, Refinement_precision, Max_iter, tol_fact)
end


end
