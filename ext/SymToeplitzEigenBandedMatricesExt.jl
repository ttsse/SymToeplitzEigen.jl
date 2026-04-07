module SymToeplitzEigenBandedMatricesExt

using LinearAlgebra
using SparseArrays
using SymToeplitzEigen
using BandedMatrices

function SymToeplitzEigen.EigenRef(
    A::BandedMatrices.AbstractBandedMatrix{T};
    banded_backend::Symbol=:auto,
    banded_density_threshold::Float64=0.35,
    banded_bandwidth_threshold::Float64=0.5,
    kwargs...
) where {T <: Union{Float32, Float64}}
    if size(A, 1) != size(A, 2)
        throw(ArgumentError("Given matrix must be square."))
    end
    if !issymmetric(A)
        throw(ArgumentError("Given banded matrix must be symmetric."))
    end
    if !(banded_backend in (:auto, :dense, :sparse))
        throw(ArgumentError("Invalid banded_backend=$(banded_backend). Expected :auto, :dense, or :sparse."))
    end
    if !(0.0 < banded_density_threshold <= 1.0)
        throw(ArgumentError("banded_density_threshold must satisfy 0 < threshold <= 1."))
    end
    if !(0.0 < banded_bandwidth_threshold <= 1.0)
        throw(ArgumentError("banded_bandwidth_threshold must satisfy 0 < threshold <= 1."))
    end

    if banded_backend == :dense
        return SymToeplitzEigen.EigenRef(Matrix(A); kwargs...)
    end

    As = sparse(A)
    sparse_kwargs = merge((sparse_bandwidth_threshold=banded_bandwidth_threshold,), kwargs)
    if banded_backend == :sparse
        return SymToeplitzEigen.EigenRef(As; sparse_kwargs...)
    end

    n = size(A, 1)
    density = nnz(As) / (n * n)
    lband, uband = BandedMatrices.bandwidths(A)
    bandwidth_ratio = (lband + uband + 1) / n
    solve_mode = get(kwargs, :solve_mode, :fast)

    if density <= banded_density_threshold || (solve_mode == :fast && bandwidth_ratio <= banded_bandwidth_threshold)
        return SymToeplitzEigen.EigenRef(As; sparse_kwargs...)
    end

    return SymToeplitzEigen.EigenRef(Matrix(A); kwargs...)
end

function SymToeplitzEigen.EigenRef(
    A::BandedMatrices.AbstractBandedMatrix{BigFloat};
    banded_backend::Symbol=:auto,
    banded_density_threshold::Float64=0.35,
    banded_bandwidth_threshold::Float64=0.5,
    kwargs...
)
    if size(A, 1) != size(A, 2)
        throw(ArgumentError("Given matrix must be square."))
    end
    if !issymmetric(A)
        throw(ArgumentError("Given banded matrix must be symmetric."))
    end
    if !(banded_backend in (:auto, :dense, :sparse))
        throw(ArgumentError("Invalid banded_backend=$(banded_backend). Expected :auto, :dense, or :sparse."))
    end
    if !(0.0 < banded_density_threshold <= 1.0)
        throw(ArgumentError("banded_density_threshold must satisfy 0 < threshold <= 1."))
    end
    if !(0.0 < banded_bandwidth_threshold <= 1.0)
        throw(ArgumentError("banded_bandwidth_threshold must satisfy 0 < threshold <= 1."))
    end

    if banded_backend == :dense
        return SymToeplitzEigen.EigenRef(Matrix{BigFloat}(A); kwargs...)
    end

    A64_sparse = Float64.(sparse(A))
    sparse_kwargs = merge((sparse_bandwidth_threshold=banded_bandwidth_threshold,), kwargs)

    if banded_backend == :sparse
        return SymToeplitzEigen.EigenRef(A64_sparse; sparse_kwargs...)
    end

    n = size(A, 1)
    density = nnz(A64_sparse) / (n * n)
    lband, uband = BandedMatrices.bandwidths(A)
    bandwidth_ratio = (lband + uband + 1) / n
    solve_mode = get(kwargs, :solve_mode, :fast)

    if density <= banded_density_threshold || (solve_mode == :fast && bandwidth_ratio <= banded_bandwidth_threshold)
        return SymToeplitzEigen.EigenRef(A64_sparse; sparse_kwargs...)
    end

    return SymToeplitzEigen.EigenRef(Matrix{BigFloat}(A); kwargs...)
end

end
