module SymToeplitzEigenToeplitzMatricesExt

using LinearAlgebra
using SymToeplitzEigen
using ToeplitzMatrices

function _first_column(A::AbstractMatrix)
    n = size(A, 1)
    c = Vector{eltype(A)}(undef, n)
    @inbounds begin
        for kk in 1:n
            c[kk] = A[kk, 1]
        end
    end
    return c
end

function SymToeplitzEigen.EigenRef(A::ToeplitzMatrices.SymmetricToeplitz{T}; kwargs...) where {T <: Real}
    c = _first_column(A)
    return SymToeplitzEigen.EigenRef(length(c), c; kwargs...)
end

function SymToeplitzEigen.EigenRef(A::ToeplitzMatrices.Toeplitz{T}; kwargs...) where {T <: Real}
    if !issymmetric(A)
        throw(ArgumentError("EigenRef requires symmetric input. Use SymmetricToeplitz or a symmetric Toeplitz matrix."))
    end

    c = _first_column(A)
    return SymToeplitzEigen.EigenRef(length(c), c; kwargs...)
end

end
