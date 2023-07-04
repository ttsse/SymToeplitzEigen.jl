
"""
    Compute_CUDA_Eigen(Tn)
    
Computes the eigenvalues (and vectors) of a symmetric matrix Tn using the CUDA implementation of syevd with Julia wrapper from CUDA.jl
"""
function Compute_CUDA_Eigen(Tn :: Array)
    Tn_cuda = CuArray(Tn)
    result = CUDA.@sync CUDA.CUSOLVER.syevd!('V','U', Tn_cuda)
    return Array(result[1]), Array(result[2])
end