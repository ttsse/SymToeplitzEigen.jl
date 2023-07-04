module SymToeplitzEigen

using LinearAlgebra

include("refinement.jl")
include("helper.jl")

export EigenRef

"""
    EigenRef(Tn; Low_pres_type = Float32, Refinement_precision = 256, Max_iter :: Integer = 10, tol_fact :: Integer = 1)

Computes the eigenvalues (and vectors) of the symmetric Toeplitz matrix `Tn`. The eigenvalues (and vectors) are first computed in "low" precision (specified using Low_pres_type with values of either Float32 or Float64). The calculated eigenvalues (and vectors) are then improved using an Iterative Refinement procedure in parallel up to a precision `Refinement_precision`, with a maximum of `Max_iter` iterations and with an error tolerance of `tol_fact*eps(BigFloat)`.

---
    EigenRef(n, vc; Low_pres_type = Float32, Refinement_precision = 256, Max_iter :: Integer = 10, tol_fact :: Integer = 1)    

Convenience function that will construct the symmetric Toeplitz matrix of size n x n with the first column given by vector vc and then return the computed eigenvalues and eigenvectors.
"""
    function EigenRef(Tn :: Array{T}; Low_pres_type :: Type = Float32, Refinement_precision :: Integer = 256, Max_iter :: Integer = 100, tol_fact :: Integer = 1) where {T <: Union{Float64, Float32}}
        return _Eigen(Tn, Refinement_precision, Max_iter, tol_fact)
    end


    function EigenRef(Tn :: Array{BigFloat}; Low_pres_type :: Type = Float32, Refinement_precision :: Integer = 256, Max_iter :: Integer = 100, tol_fact :: Integer = 1)
        Tn_low = Low_pres_type.(Tn)
        return _Eigen(Tn_low, Tn[1].prec, Max_iter, tol_fact)
    end


    function EigenRef(n :: Integer, vc :: Array; Low_pres_type :: Type = Float32, Refinement_precision :: Integer = 256, Max_iter :: Integer = 100, tol_fact :: Integer = 1)
        Tn = toeplitz(n, Low_pres_type.(vc), Low_pres_type.(vc))
        return _Eigen(Tn, Refinement_precision, Max_iter, tol_fact)
    end

    function _Eigen(Tn :: Array{T}, Refinement_precision :: Integer, Max_iter :: Integer, tol_fact :: Integer) where {T <: Union{Float64, Float32}}
        vals, vecs = eigen(Tn)
        return Refinement(Tn, vals, vecs, Refinement_precision, Max_iter, tol_fact)
    end


end #SymToeplitzEigen
