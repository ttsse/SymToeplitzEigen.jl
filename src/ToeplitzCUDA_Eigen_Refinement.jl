module ToeplitzCUDA_Eigen_Refinement

using CUDA

using IterativeRefinement

using LinearAlgebra

export compute_Toeplitz_Eigen

"""
    compute_Toeplitz_Eigen(n, vc; CUDA_type = Float32, Refinement_precision = 256, Refinement_Iterations = :single)

Computes the eigenvalues (and vectors) of the symmetric Toeplitz matrix of size n x n with the first column given by vector vc. The eigenvalues (and vectors) are first computed in "low" precision (specified using CUDA_type with values of either Float32 or Float64) on the GPU. The calculated values are then improved using an Iterative Refinement procedure in parallel.
"""
function compute_Toeplitz_Eigen(n :: Integer, vc :: Array; CUDA_type :: Type = Float32, Refinement_precision :: Integer = 256, Multiple_refinement :: Bool = false, Max_iter = 5)

    Tn = toeplitz(n, CUDA_type.(vc), CUDA_type.(vc))

    CUDAvals, CUDAvecs = Compute_CUDA_Eigen(Tn)

    return Refinement(Tn, CUDAvals, CUDAvecs, Refinement_precision, Multiple_refinement, Max_iter)

end

"""
    toeplitz(n, vc, vr)

Constructs a Toeplitz matrix of size n x n with the first column given by vector vc and first row given by vr. The element type of the given vectors determine the element type of the Toeplitz matrix.
"""
function toeplitz(n :: Integer, vc :: Array{T, 1}, vr :: Array{T, 1}) where T <: Number
    Tn = zeros(T,n,n)
    fill!.(view.((Tn,), [diagind(Tn, 1-ii) for ii in eachindex(vc)]), view(vc, :))
    fill!.(view.((Tn,), [diagind(Tn, jj-1) for jj in Iterators.drop(eachindex(vr), 1)]), @view vr[2:end])
    return Tn
end

# function toeplitzBanded(n :: Integer, vc :: Array{T, 1}, vr :: Array{T, 1}) where T <: Number
#     Tn = BandedMatrix(Zeros(T,n,n), (length(vc)-1, length(vr)-1))
#     fill!.(view.((Tn,), [diagind(Tn, 1-ii) for ii in eachindex(vc)]), view(vc, :))
#     fill!.(view.((Tn,), [diagind(Tn, jj-1) for jj in Iterators.drop(eachindex(vr), 1)]), @view vr[2:end])
#     return Tn
# end

function Compute_CUDA_Eigen(Tn :: Array)
    Tn_cuda = CuArray(Tn)
    result = CUDA.@sync CUDA.CUSOLVER.syevd!('V','U', Tn_cuda)
    return Array(result[1]), Array(result[2])
end

function Refinement(Tn :: Array, CUDAvals :: Array, CUDAvecs :: Array, Refinement_precision :: Integer, Multiple_refinement :: Bool, Max_iter :: Integer)
    setprecision(BigFloat, Refinement_precision)
    n = length(CUDAvals)
    refinedVals = zeros(BigFloat,n)
    refinedVecs = zeros(BigFloat,n,n)
    Threads.@threads for ii = 1:n
        newVal, newVec = rfeigen(Tn, BigFloat.(view(CUDAvecs, :, ii)), BigFloat(CUDAvals[ii]), BigFloat)
        if Multiple_refinement
            err_norm = norm((Tn - newVal*I)*newVec)
            iter = 0
            while err_norm > 100*eps(BigFloat) && iter < Max_iter
                iter += 1
                newVal, newVec = rfeigen(Tn, newVec, newVal, BigFloat)
                err_norm = norm((Tn - newVal*I)*newVec)
            end
        end
        refinedVals[ii] = newVal
        refinedVecs[:,ii] .= newVec
    end
    return refinedVals, refinedVecs
end

end
