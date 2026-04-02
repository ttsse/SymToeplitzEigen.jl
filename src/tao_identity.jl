"""
    _tao_safe_gap(vals, i, gap_factor)

Returns `(is_safe, min_gap, gap_tol)` for Tao-identity use on pair `i`.
"""
function _tao_safe_gap(vals::AbstractVector{BigFloat}, i::Integer, gap_factor::Float64)
    λ = vals[i]
    min_gap = BigFloat(Inf)
    @inbounds for k in eachindex(vals)
        if k != i
            d = abs(λ - vals[k])
            if d < min_gap
                min_gap = d
            end
        end
    end

    gap_tol = BigFloat(gap_factor) * eps(BigFloat) * max(abs(λ), one(BigFloat))
    return min_gap > gap_tol, min_gap, gap_tol
end

"""
    tao_component_magnitude_squared(A, vals, i, j; gap_factor=1e4)

Computes Tao's identity estimate of `|v[i,j]|^2` for symmetric matrices.
Returns `nothing` when safety gates fail.
"""
function tao_component_magnitude_squared(
    A::AbstractMatrix{T},
    vals::AbstractVector,
    i::Integer,
    j::Integer;
    gap_factor::Float64=1e4,
) where {T <: Real}
    n = size(A, 1)
    if size(A, 1) != size(A, 2)
        throw(ArgumentError("A must be square"))
    end
    if !issymmetric(A)
        return nothing
    end
    if !(1 <= i <= length(vals)) || !(1 <= j <= n)
        throw(ArgumentError("Indices out of bounds"))
    end

    vals_big = BigFloat.(real.(vals))
    safe, _, _ = _tao_safe_gap(vals_big, i, gap_factor)
    if !safe
        return nothing
    end

    λ = vals_big[i]
    idx = collect(1:n)
    deleteat!(idx, j)

    # LinearAlgebra does not provide a stable eigvals path for generic BigFloat dense matrices,
    # so use a symmetric Float64 solve for the minor and lift back to BigFloat for products.
    A_minor = Matrix{Float64}(A[idx, idx])
    μ = eigvals(Symmetric(A_minor))

    numer = one(BigFloat)
    for u in μ
        numer *= (λ - BigFloat(real(u)))
    end

    denom = one(BigFloat)
    for k in eachindex(vals_big)
        if k != i
            denom *= (λ - vals_big[k])
        end
    end

    if abs(denom) <= eps(BigFloat)
        return nothing
    end

    return abs(numer / denom)
end

"""
    TaoIdentityReport(A, vals, vecs; gap_factor=1e4, max_pairs=10, max_components=3)

Post-refinement diagnostic using Tao's identity for symmetric matrices.
"""
function TaoIdentityReport(
    A::AbstractMatrix{T},
    vals::AbstractVector,
    vecs::AbstractMatrix;
    gap_factor::Float64=1e4,
    max_pairs::Integer=10,
    max_components::Integer=3,
) where {T <: Real}
    n = size(A, 1)
    if size(A, 1) != size(A, 2)
        return (enabled=false, reason="matrix-not-square")
    end
    if !issymmetric(A)
        return (enabled=false, reason="matrix-not-symmetric")
    end

    vals_big = BigFloat.(real.(vals))
    vecs_big = BigFloat.(vecs)

    pairs = min(length(vals_big), max_pairs)
    checked_pairs = 0
    skipped_pairs = 0
    checked_components = 0

    relerrs = BigFloat[]
    details = NamedTuple[]

    for i in 1:pairs
        safe, min_gap, gap_tol = _tao_safe_gap(vals_big, i, gap_factor)
        if !safe
            skipped_pairs += 1
            push!(details, (pair=i, skipped=true, reason="small-gap", min_gap=min_gap, gap_tol=gap_tol))
            continue
        end

        checked_pairs += 1
        cands = unique([argmax(abs.(view(vecs_big, :, i))), 1, n])
        comps = cands[1:min(length(cands), max_components)]

        for j in comps
            pred = tao_component_magnitude_squared(A, vals_big, i, j; gap_factor=gap_factor)
            if isnothing(pred)
                continue
            end
            actual = abs(vecs_big[j, i])^2
            rel = abs(pred - actual) / max(actual, eps(BigFloat))
            checked_components += 1
            push!(relerrs, rel)
            push!(details, (pair=i, component=j, predicted=pred, actual=actual, relative_error=rel, skipped=false))
        end
    end

    mean_rel = isempty(relerrs) ? BigFloat(Inf) : sum(relerrs) / length(relerrs)
    max_rel = isempty(relerrs) ? BigFloat(Inf) : maximum(relerrs)

    return (
        enabled=true,
        reason="ok",
        checked_pairs=checked_pairs,
        skipped_pairs=skipped_pairs,
        checked_components=checked_components,
        mean_relative_error=mean_rel,
        max_relative_error=max_rel,
        details=details,
    )
end

"""
    TaoScaleEigenvectors!(A, vals, vecs; gap_factor=1e4, max_pairs=0)

Experimental initialization/scaling prototype: rescales each selected eigenvector
so its dominant component magnitude aligns with Tao's identity estimate.
"""
function TaoScaleEigenvectors!(
    A::AbstractMatrix{T},
    vals::AbstractVector,
    vecs::AbstractMatrix;
    gap_factor::Float64=1e4,
    max_pairs::Integer=0,
) where {T <: Real}
    if !issymmetric(A)
        return (scaled_pairs=0, skipped_pairs=length(vals), reason="matrix-not-symmetric")
    end

    m = length(vals)
    p = max_pairs > 0 ? min(max_pairs, m) : m

    scaled_pairs = 0
    skipped_pairs = 0

    for i in 1:p
        j = argmax(abs.(view(vecs, :, i)))
        pred = tao_component_magnitude_squared(A, vals, i, j; gap_factor=gap_factor)
        if isnothing(pred)
            skipped_pairs += 1
            continue
        end

        actual = abs(BigFloat(vecs[j, i]))^2
        if actual <= eps(BigFloat)
            skipped_pairs += 1
            continue
        end

        scale = sqrt(pred / actual)
        vecs[:, i] .*= scale
        scaled_pairs += 1
    end

    return (scaled_pairs=scaled_pairs, skipped_pairs=skipped_pairs, reason="ok")
end
