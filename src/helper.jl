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

const _BIGFLOAT_MATVEC_BLOCK = 32
const _BIGFLOAT_FILL_BLOCK = 32

"""
    normalize_Infnorm!(x)

Normalises the input array `x` wrt the maximum norm (normalised so that maximal value in the array is 1) and returns the index of the maximal value.
"""
@inline function normalize_Infnorm!(x :: AbstractVector)
    n = length(x)
    s=1
    xm = abs(x[1])
    for j in 1:n
        t = abs(x[j])
        if t > xm
            xm = t
            s = j
        end
    end
    # my_div_vec!(n, x, xm)
    x ./= xm
    return s
end

"""
    y_div_vec!(n, V, val)

Divides each element of vector `V` of length `n` with the value `val`.
"""
@inline function my_div_vec!(n::Integer, V :: StridedVector{BigFloat}, val :: BigFloat)
    @inbounds begin
        for kk in 1:n
            my_div!(V[kk], V[kk], val)
        end
    end
end

"""
    my_mat_vec_mul!(n, R, A, B)

Computes the matrix vector multiplication `R = -A*B` where `A` is a matrix and `B` is a vector.
"""
@inline function my_neg_mat_vec_mul!(n::Integer, R :: StridedVector{BigFloat}, A :: AbstractArray{BigFloat}, B :: StridedVector{BigFloat})
    @inbounds begin
        for kk in 1:n
            setzero!(R[kk])
            for j0 in 1:_BIGFLOAT_MATVEC_BLOCK:n
                j1 = min(j0 + _BIGFLOAT_MATVEC_BLOCK - 1, n)
                for jj in j0:j1
                    my_fma!(R[kk], R[kk], A[kk, jj], B[jj])
                end
            end
            neg!(R[kk])
        end
    end
end

"""
    my_neg_sparse_mat_vec_mul!(n, R, A, x, x64, y64)

Computes `R = -A*x` for sparse `A` while reusing Float64 workspaces.
The inputs `x64` and `y64` are thread-local scratch vectors.
"""
@inline function my_neg_sparse_mat_vec_mul!(
    n::Integer,
    R :: StridedVector{BigFloat},
    A :: SparseMatrixCSC{Float64, Ti},
    x :: StridedVector{BigFloat},
    x64 :: StridedVector{Float64},
    y64 :: StridedVector{Float64},
) where {Ti <: Integer}
    my_vec_demote!(n, x64, x)
    mul!(y64, A, x64)

    @inbounds begin
        for kk in 1:n
            setvalue_F!(R[kk], -y64[kk])
        end
    end
end

"""
    my_neg_sym_toeplitz_mat_vec!(n, R, c, x)

Computes `R = -T*x` for a symmetric Toeplitz matrix `T` represented by
its first column `c`. If `length(c) < n`, trailing diagonals are treated as zero.
"""
@inline function my_neg_sym_toeplitz_mat_vec!(n::Integer, R :: StridedVector{BigFloat}, c :: StridedVector{BigFloat}, x :: StridedVector{BigFloat})
    m = length(c)
    @inbounds begin
        if m >= n
            for kk in 1:n
                setzero!(R[kk])
                for j0 in 1:_BIGFLOAT_MATVEC_BLOCK:n
                    j1 = min(j0 + _BIGFLOAT_MATVEC_BLOCK - 1, n)
                    for jj in j0:j1
                        my_fma!(R[kk], R[kk], c[abs(kk - jj) + 1], x[jj])
                    end
                end
                neg!(R[kk])
            end
        else
            for kk in 1:n
                setzero!(R[kk])
                for j0 in 1:_BIGFLOAT_MATVEC_BLOCK:n
                    j1 = min(j0 + _BIGFLOAT_MATVEC_BLOCK - 1, n)
                    for jj in j0:j1
                        idx = abs(kk - jj) + 1
                        if idx <= m
                            my_fma!(R[kk], R[kk], c[idx], x[jj])
                        end
                    end
                end
                neg!(R[kk])
            end
        end
    end
end

"""
    my_neg_sym_toeplitz_mat_vec_fft!(n, R, c, x, x64, fft_cache)

Optional FFTW-backed Toeplitz mat-vec kernel for `R = -T*x`, where `c` is the
Float64 first column of `T` and `x` is BigFloat. This method is provided by the
FFTW extension when FFTW.jl is loaded.
"""
@inline function my_neg_sym_toeplitz_mat_vec_fft!(
    n::Integer,
    R :: AbstractVector{BigFloat},
    c :: AbstractVector{<:Real},
    x :: AbstractVector{BigFloat},
    x64 :: AbstractVector{Float64},
    fft_cache :: Base.RefValue,
)
    throw(ArgumentError("toeplitz_kernel=:fft requires FFTW.jl to be loaded. Add FFTW.jl and run `using FFTW` before calling EigenRef."))
end

"""
    my_fill_sym_toeplitz!(n, A, c)

Fills matrix `A` as an `n x n` symmetric Toeplitz matrix from first column `c`.
If `length(c) < n`, trailing diagonals are filled with zeros.
"""
@inline function my_fill_sym_toeplitz!(n::Integer, A :: AbstractMatrix{Float64}, c :: StridedVector{Float64})
    m = length(c)
    @inbounds begin
        if m >= n
            for j0 in 1:_BIGFLOAT_FILL_BLOCK:n
                j1 = min(j0 + _BIGFLOAT_FILL_BLOCK - 1, n)
                for jj in j0:j1
                    for kk in 1:n
                        A[kk, jj] = c[abs(kk - jj) + 1]
                    end
                end
            end
        else
            for j0 in 1:_BIGFLOAT_FILL_BLOCK:n
                j1 = min(j0 + _BIGFLOAT_FILL_BLOCK - 1, n)
                for jj in j0:j1
                    for kk in 1:n
                        idx = abs(kk - jj) + 1
                        A[kk, jj] = idx <= m ? c[idx] : 0.0
                    end
                end
            end
        end
    end
end

@inline function my_fill_sym_toeplitz!(n::Integer, A :: AbstractMatrix{BigFloat}, c :: StridedVector{BigFloat})
    m = length(c)
    @inbounds begin
        if m >= n
            for j0 in 1:_BIGFLOAT_FILL_BLOCK:n
                j1 = min(j0 + _BIGFLOAT_FILL_BLOCK - 1, n)
                for jj in j0:j1
                    for kk in 1:n
                        setvalue_BF!(A[kk, jj], c[abs(kk - jj) + 1])
                    end
                end
            end
        else
            for j0 in 1:_BIGFLOAT_FILL_BLOCK:n
                j1 = min(j0 + _BIGFLOAT_FILL_BLOCK - 1, n)
                for jj in j0:j1
                    for kk in 1:n
                        idx = abs(kk - jj) + 1
                        if idx <= m
                            setvalue_BF!(A[kk, jj], c[idx])
                        else
                            setzero!(A[kk, jj])
                        end
                    end
                end
            end
        end
    end
end

"""
    my_add_diag_elements!(n, A, val)

Adds the value `val` to the diagonal elements of the matrix `A` of length `n`.
"""
@inline function my_add_diag_elements!(n::Integer, A :: AbstractArray{BigFloat}, val :: BigFloat)
    @inbounds begin
        for kk in 1:n
            my_add!(A[kk, kk], A[kk, kk], val)
        end
    end
end

"""
    my_vec_setvalue_prom!(n, A, B)

Sets the vector `A` equal to the elements of vector `B` where we promote the type from `Float64` to `BigFloat`.
"""
@inline function my_vec_setvalue_prom!(n::Integer, A :: StridedVector{BigFloat}, B :: StridedVector{T}) where {T <: Union{Float32, Float64}}
    @inbounds begin
        for kk in 1:n
            setvalue_F!(A[kk], B[kk])
        end
    end
end

"""
my_vec_setvalue!(n, A, B)

Sets the vector `A` equal to the elements of vector `B`.
"""
@inline function my_vec_setvalue!(n::Integer, A :: StridedVector{BigFloat}, B :: StridedVector{BigFloat})
    @inbounds begin
        for kk in 1:n
            setvalue_BF!(A[kk], B[kk])
        end
    end
end

"""
    my_mat_setvalue!(n, A, B)

Sets matrix `A` equal to matrix `B` for square `n x n` BigFloat matrices using in-place writes.
"""
@inline function my_mat_setvalue!(n::Integer, A :: AbstractMatrix{BigFloat}, B :: AbstractMatrix{BigFloat})
    @inbounds begin
        for j0 in 1:_BIGFLOAT_FILL_BLOCK:n
            j1 = min(j0 + _BIGFLOAT_FILL_BLOCK - 1, n)
            for jj in j0:j1
                for kk in 1:n
                    setvalue_BF!(A[kk, jj], B[kk, jj])
                end
            end
        end
    end
end

"""
    my_update_shifted_pivot_system!(n, B, lambda_prev, lambda_now, s, x)

Updates a previously factorized linear system matrix in place when only the
diagonal shift and replacement pivot column change.
"""
@inline function my_update_shifted_pivot_system!(
    n::Integer,
    B :: AbstractMatrix{Float64},
    lambda_prev :: Float64,
    lambda_now :: Float64,
    s :: Integer,
    x :: StridedVector{Float64},
)
    delta = lambda_prev - lambda_now
    @inbounds begin
        for kk in 1:n
            B[kk, kk] += delta
            B[kk, s] = -x[kk]
        end
    end
end

@inline function my_update_shifted_pivot_system!(
    n::Integer,
    B :: AbstractMatrix{BigFloat},
    lambda_prev :: BigFloat,
    lambda_now :: BigFloat,
    s :: Integer,
    x :: StridedVector{BigFloat},
)
    delta = lambda_prev - lambda_now
    @inbounds begin
        for kk in 1:n
            my_add!(B[kk, kk], B[kk, kk], delta)
            setvalue_BF!(B[kk, s], x[kk])
            neg!(B[kk, s])
        end
    end
end

"""
    my_fma!(w, z, x, y)

Non-allocating fma (`w = z + x*y`) for `BigFloat`.
"""
function my_fma!(w::BigFloat, z::BigFloat, x::BigFloat, y::BigFloat)
    ccall(("mpfr_fma",Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), w, x, y, z, Base.MPFR.MPFRRoundNearest)
end

"""
    my_add!(z, y, x)

Non-allocating add (`z = x + y`) for `BigFloat`.
"""
function my_add!(z::BigFloat, x::BigFloat, y::BigFloat)
    ccall(("mpfr_add",Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, y, Base.MPFR.MPFRRoundNearest)
end

"""
    my_div!(z, y, x)

Non-allocating div (`z = x / y`) for `BigFloat`.
"""
function my_div!(z::BigFloat, x::BigFloat, y::BigFloat)
    ccall(("mpfr_div",Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, y, Base.MPFR.MPFRRoundNearest)
end

"""
    neg!(z)

Non-allocating sets `z = -z`.
"""
function neg!(z::BigFloat)
    ccall(("mpfr_neg", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, z, Base.MPFR.MPFRRoundNearest)
end

"""
    setzero!(z)

Sets the value of `z` equal to 0.
"""
function setzero!(z::BigFloat)
    ccall(("mpfr_set_ui", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Culong, Base.MPFR.MPFRRoundingMode), z, 0, Base.MPFR.MPFRRoundNearest)
end

"""
    setvalue_BF!(z, x)

Sets the value of `z` equal to `x` with both `BigFloat`.
"""
function setvalue_BF!(z::BigFloat, x::BigFloat)
    ccall(("mpfr_set", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, Base.MPFR.MPFRRoundNearest)
end

"""
    setvalue_F!(z, x)

Sets the value of `z` equal to `x` with `x` being promoted from `Float32` (or `Float64`) to `BigFloat`.
"""
function setvalue_F!(z::BigFloat, x::Float32)
    ccall(("mpfr_set_d", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Cfloat, Base.MPFR.MPFRRoundingMode), z, x, Base.MPFR.MPFRRoundNearest)
end

function setvalue_F!(z::BigFloat, x::Float64)
    ccall(("mpfr_set_d", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Cdouble, Base.MPFR.MPFRRoundingMode), z, x, Base.MPFR.MPFRRoundNearest)
end

"""
    my_add_vec!(n, A, B)

Sets `A .= A .+ B` for vectors of length `n`.
"""
@inline function my_add_vec!(n::Integer, A :: StridedVector{BigFloat}, B :: StridedVector{BigFloat})
    @inbounds begin
        for kk in 1:n
            my_add!(A[kk], A[kk], B[kk])
        end
    end
end

"""
    my_add_scaled_vec!(n, A, scale, B)

Sets `A .= A .+ scale .* B` for vectors of length `n`.
"""
@inline function my_add_scaled_vec!(n::Integer, A :: StridedVector{BigFloat}, scale :: BigFloat, B :: StridedVector{BigFloat})
    @inbounds begin
        for kk in 1:n
            my_fma!(A[kk], A[kk], scale, B[kk])
        end
    end
end

"""
    my_vec_demote!(n, A, B)

Sets `A` equal to `B` with each element of `B` converted from `BigFloat` to `Float64`.
"""
@inline function my_vec_demote!(n::Integer, A :: StridedVector{Float64}, B :: StridedVector{BigFloat})
    @inbounds begin
        for kk in 1:n
            A[kk] = Float64(B[kk])
        end
    end
end

"""
    my_maxabs(n, x)

Returns the infinity norm of vector `x` (maximum absolute value).
"""
@inline function my_maxabs(n::Integer, x :: StridedVector{BigFloat})
    xm = abs(x[1])
    @inbounds begin
        for kk in 2:n
            t = abs(x[kk])
            if t > xm
                xm = t
            end
        end
    end
    return xm
end
