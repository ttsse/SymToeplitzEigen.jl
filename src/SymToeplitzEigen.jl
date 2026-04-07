module SymToeplitzEigen

using LinearAlgebra
using OhMyThreads
using SparseArrays

include("refinement.jl")
include("helper.jl")
include("sparse_refinement.jl")
include("progress.jl")
include("nonsymmetric.jl")

export EigenRef, EigenRefSparse, EigenRefNonSym

mutable struct _ToeplitzCacheEntry{T}
    n::Int
    vc::Vector{T}
    mat::Matrix{T}
end

const _toeplitz_cache_ref = Ref{Any}(nothing)
const _toeplitz_cache_lock = ReentrantLock()

function _cached_toeplitz(n::Integer, vc::Vector{T}; reuse::Bool=true) where T <: Number
    if !reuse
        return toeplitz(n, vc, vc)
    end

    lock(_toeplitz_cache_lock)
    try
        entry = _toeplitz_cache_ref[]
        if entry isa _ToeplitzCacheEntry{T}
            if entry.n == n && entry.vc == vc
                return entry.mat
            end
        end

        Tn = toeplitz(n, vc, vc)
        _toeplitz_cache_ref[] = _ToeplitzCacheEntry{T}(n, copy(vc), Tn)
        return Tn
    finally
        unlock(_toeplitz_cache_lock)
    end
end

"""
    EigenRef(Tn; Low_pres_type = Float64, Refinement_precision = 256, Max_iter :: Integer = 10, tol_fact :: Integer = 1, solve_mode = :fast, stall_ratio = 0.9, stall_iters = 3, return_status = false, adaptive_precision_escalation = false, escalation_precision = 0, escalation_extra_iter = 10, toeplitz_kernel = :auto, toeplitz_auto_threshold = 96, show_progress = false)

Computes the eigenvalues (and vectors) of the symmetric Toeplitz matrix `Tn`. The eigenvalues (and vectors) are first computed in "low" precision (specified using Low_pres_type with values of either Float32 or Float64). The calculated eigenvalues (and vectors) are then improved using an Iterative Refinement procedure in parallel up to a precision `Refinement_precision`, with a maximum of `Max_iter` iterations and with an error tolerance of `tol_fact*eps(BigFloat)`.

---
    EigenRef(n, vc; Low_pres_type = Float64, Refinement_precision = 256, Max_iter :: Integer = 10, tol_fact :: Integer = 1, solve_mode = :fast, stall_ratio = 0.9, stall_iters = 3, return_status = false, adaptive_precision_escalation = false, escalation_precision = 0, escalation_extra_iter = 10, toeplitz_kernel = :auto, toeplitz_auto_threshold = 96, reuse_toeplitz_cache = true, show_progress = false)

Convenience function that will construct the symmetric Toeplitz matrix of size n x n with the first column given by vector vc and then return the computed eigenvalues and eigenvectors.
"""
    function EigenRef(Tn :: Array{T}; Low_pres_type :: Type = Float64, Refinement_precision :: Integer = 256, Max_iter :: Integer = 100, tol_fact :: Integer = 1, solve_mode :: Symbol = :fast, stall_ratio :: Float64 = 0.9, stall_iters :: Integer = 3, return_status :: Bool = false, adaptive_precision_escalation :: Bool = false, escalation_precision :: Integer = 0, escalation_extra_iter :: Integer = 10, toeplitz_kernel :: Symbol = :auto, toeplitz_auto_threshold :: Integer = 96, show_progress :: Bool = false) where {T <: Union{Float64, Float32}}
        @assert issymmetric(Tn) "Given matrix must be symmetric"
        return _Eigen(Tn, Refinement_precision, Max_iter, tol_fact; solve_mode=solve_mode, stall_ratio=stall_ratio, stall_iters=stall_iters, return_status=return_status, adaptive_precision_escalation=adaptive_precision_escalation, escalation_precision=escalation_precision, escalation_extra_iter=escalation_extra_iter, toeplitz_kernel=toeplitz_kernel, toeplitz_auto_threshold=toeplitz_auto_threshold, show_progress=show_progress)
    end


    function EigenRef(Tn :: Array{BigFloat}; Low_pres_type :: Type = Float64, Refinement_precision :: Integer = 256, Max_iter :: Integer = 100, tol_fact :: Integer = 1, solve_mode :: Symbol = :fast, stall_ratio :: Float64 = 0.9, stall_iters :: Integer = 3, return_status :: Bool = false, adaptive_precision_escalation :: Bool = false, escalation_precision :: Integer = 0, escalation_extra_iter :: Integer = 10, toeplitz_kernel :: Symbol = :auto, toeplitz_auto_threshold :: Integer = 96, show_progress :: Bool = false)
        Tn_low = Low_pres_type.(Tn)
        return _Eigen(Tn_low, Tn[1].prec, Max_iter, tol_fact; solve_mode=solve_mode, stall_ratio=stall_ratio, stall_iters=stall_iters, return_status=return_status, adaptive_precision_escalation=adaptive_precision_escalation, escalation_precision=escalation_precision, escalation_extra_iter=escalation_extra_iter, toeplitz_kernel=toeplitz_kernel, toeplitz_auto_threshold=toeplitz_auto_threshold, show_progress=show_progress)
    end


    function EigenRef(n :: Integer, vc :: Array; Low_pres_type :: Type = Float64, Refinement_precision :: Integer = 256, Max_iter :: Integer = 100, tol_fact :: Integer = 1, solve_mode :: Symbol = :fast, stall_ratio :: Float64 = 0.9, stall_iters :: Integer = 3, return_status :: Bool = false, adaptive_precision_escalation :: Bool = false, escalation_precision :: Integer = 0, escalation_extra_iter :: Integer = 10, toeplitz_kernel :: Symbol = :auto, toeplitz_auto_threshold :: Integer = 96, reuse_toeplitz_cache :: Bool = true, show_progress :: Bool = false)
        vc_low = Low_pres_type.(vc)
        Tn = _cached_toeplitz(n, vc_low; reuse=reuse_toeplitz_cache)
        return _Eigen(Tn, Refinement_precision, Max_iter, tol_fact; solve_mode=solve_mode, stall_ratio=stall_ratio, stall_iters=stall_iters, return_status=return_status, adaptive_precision_escalation=adaptive_precision_escalation, escalation_precision=escalation_precision, escalation_extra_iter=escalation_extra_iter, toeplitz_kernel=toeplitz_kernel, toeplitz_column=vc_low, toeplitz_auto_threshold=toeplitz_auto_threshold, show_progress=show_progress)
    end

    function _Eigen(Tn :: Array{T}, Refinement_precision :: Integer, Max_iter :: Integer, tol_fact :: Integer; solve_mode :: Symbol = :fast, stall_ratio :: Float64 = 0.9, stall_iters :: Integer = 3, return_status :: Bool = false, adaptive_precision_escalation :: Bool = false, escalation_precision :: Integer = 0, escalation_extra_iter :: Integer = 10, toeplitz_kernel :: Symbol = :auto, toeplitz_column = nothing, toeplitz_auto_threshold :: Integer = 96, show_progress :: Bool = false) where {T <: Union{Float64, Float32}}

        vals, vecs = eigen(Tn)

        if return_status
            return Refinement(Tn, vals, vecs, Refinement_precision, Max_iter, tol_fact; solve_mode=solve_mode, stall_ratio=stall_ratio, stall_iters=stall_iters, return_status=true, adaptive_precision_escalation=adaptive_precision_escalation, escalation_precision=escalation_precision, escalation_extra_iter=escalation_extra_iter, toeplitz_kernel=toeplitz_kernel, toeplitz_column=toeplitz_column, toeplitz_auto_threshold=toeplitz_auto_threshold, show_progress=show_progress)
        end

        return Refinement(Tn, vals, vecs, Refinement_precision, Max_iter, tol_fact; solve_mode=solve_mode, stall_ratio=stall_ratio, stall_iters=stall_iters, return_status=false, adaptive_precision_escalation=adaptive_precision_escalation, escalation_precision=escalation_precision, escalation_extra_iter=escalation_extra_iter, toeplitz_kernel=toeplitz_kernel, toeplitz_column=toeplitz_column, toeplitz_auto_threshold=toeplitz_auto_threshold, show_progress=show_progress)
    end


end #SymToeplitzEigen
