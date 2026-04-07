module SymToeplitzEigenFFTWExt

using FFTW
using SymToeplitzEigen

mutable struct _ToeplitzFFTCache
    n::Int
    spectrum::Vector{ComplexF64}
    xbuf::Vector{ComplexF64}
    ybuf::Vector{ComplexF64}
end

function _build_cache(c::StridedVector{Float64}, n::Int)
    m = min(length(c), n)
    embed = zeros(ComplexF64, 2 * n)

    @inbounds begin
        for kk in 1:m
            embed[kk] = ComplexF64(c[kk], 0.0)
        end

        for kk in 2:n
            embed[2 * n - kk + 2] = ComplexF64(kk <= m ? c[kk] : 0.0, 0.0)
        end
    end

    spectrum = FFTW.fft(embed)
    xbuf = zeros(ComplexF64, 2 * n)
    ybuf = zeros(ComplexF64, 2 * n)
    return _ToeplitzFFTCache(n, spectrum, xbuf, ybuf)
end

function SymToeplitzEigen.my_neg_sym_toeplitz_mat_vec_fft!(
    n::Integer,
    R::StridedVector{BigFloat},
    c::StridedVector{Float64},
    x::StridedVector{BigFloat},
    x64::StridedVector{Float64},
    fft_cache::Base.RefValue{Any},
)
    if length(c) < 1
        throw(ArgumentError("Toeplitz first-column vector must be non-empty for toeplitz_kernel=:fft."))
    end

    cache = fft_cache[]
    if !(cache isa _ToeplitzFFTCache) || cache.n != n
        cache = _build_cache(c, n)
        fft_cache[] = cache
    end
    ccache = cache::_ToeplitzFFTCache

    SymToeplitzEigen.my_vec_demote!(n, x64, x)

    fill!(ccache.xbuf, 0.0 + 0.0im)
    @inbounds begin
        for kk in 1:n
            ccache.xbuf[kk] = ComplexF64(x64[kk], 0.0)
        end
    end

    copyto!(ccache.ybuf, ccache.xbuf)
    FFTW.fft!(ccache.ybuf)
    @inbounds begin
        for kk in eachindex(ccache.ybuf)
            ccache.ybuf[kk] *= ccache.spectrum[kk]
        end
    end
    FFTW.ifft!(ccache.ybuf)

    @inbounds begin
        for kk in 1:n
            SymToeplitzEigen.setvalue_F!(R[kk], -real(ccache.ybuf[kk]))
        end
    end

    return nothing
end

end
