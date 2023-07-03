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

"""
    normalize_Infnorm!(x)

Normalises the input array `x` wrt the maximum norm (normalised so that maximal value in the array is 1) and returns the index of the maximal value.
"""
@inline function normalize_Infnorm!(x :: Array)
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
    x ./ xm
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
@inline function my_neg_mat_vec_mul!(n::Integer, R :: StridedVector{BigFloat}, A :: StridedMatrix{BigFloat}, B :: StridedVector{BigFloat})
    @inbounds begin
        for kk in 1:n
            setzero!(R[kk])
            for jj in 1:n
                my_fma!(R[kk], R[kk], A[kk, jj], B[jj]) 
            end
            neg!(R[kk])
        end
    end
end

"""
    my_add_diag_elements!(n, A, val)

Adds the value `val` to the diagonal elements of the matrix `A` of length `n`.
"""
@inline function my_add_diag_elements!(n::Integer, A :: StridedMatrix{BigFloat}, val :: BigFloat)
    @inbounds begin
        for kk in 1:n
            my_add!(A[kk, kk], A[kk, kk], val)
        end
    end
end

"""
    my_vec_setvalue_prom!(n, A, B)

Sets the vector `A` equal to the elements of vector `B` where we promote the type from Float64 to BigFloat.
"""
@inline function my_vec_setvalue_prom!(n::Integer, A :: StridedVector{BigFloat}, B :: StridedVector{Float64})
    @inbounds begin
        for kk in 1:n
            setvalue_F64!(A[kk], B[kk])
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
    my_fma!(w, z, x, y)

Non-allocating fma (`w = z + x*y`) for BigFloats.
"""
function my_fma!(w::BigFloat, z::BigFloat, x::BigFloat, y::BigFloat)
    ccall(("mpfr_fma",Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), w, x, y, z, Base.MPFR.MPFRRoundNearest)
end

"""
    my_add!(z, y, x)

Non-allocating add (`z = x + y`) for BigFloats.
"""
function my_add!(z::BigFloat, x::BigFloat, y::BigFloat)
    ccall(("mpfr_add",Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, y, Base.MPFR.MPFRRoundNearest)
end

"""
    my_div!(z, y, x)

Non-allocating div (`z = x / y`) for BigFloats.
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

Sets the value of `z` equal to 0
"""
function setzero!(z::BigFloat)
    ccall(("mpfr_set_ui", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Culong, Base.MPFR.MPFRRoundingMode), z, 0, Base.MPFR.MPFRRoundNearest)
end

"""
    setvalue_BF!(z, x)

Sets the value of `z` equal to `x` with both BigFloat
"""
function setvalue_BF!(z::BigFloat, x::BigFloat)
    ccall(("mpfr_set", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, Base.MPFR.MPFRRoundNearest)
end

"""
    setvalue_F64!(z, x)

Sets the value of `z` equal to `x` with `x` being promoted from Float64 to BigFloat
"""
function setvalue_F64!(z::BigFloat, x::Float64)
    ccall(("mpfr_set_d", Base.MPFR.libmpfr), Int32, (Ref{BigFloat}, Cdouble, Base.MPFR.MPFRRoundingMode), z, x, Base.MPFR.MPFRRoundNearest)
end