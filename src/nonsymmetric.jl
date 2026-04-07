"""
    _match_left_vectors(vals_right, vals_left, left_vectors)

Matches right-eigenvalue ordering to left-eigenvectors (computed from adjoint matrix)
using nearest conjugate eigenvalue matching.
"""
function _match_left_vectors(vals_right::AbstractVector{Complex{T}}, vals_left::AbstractVector{Complex{T}}, left_vectors::AbstractMatrix{Complex{T}}) where T
    n = length(vals_right)
    used = falses(length(vals_left))
    out = Matrix{Complex{T}}(undef, size(left_vectors, 1), n)

    @inbounds for ii in 1:n
        target = conj(vals_right[ii])
        best_idx = 1
        best_dist = Inf
        for jj in 1:length(vals_left)
            if !used[jj]
                d = abs(vals_left[jj] - target)
                if d < best_dist
                    best_dist = d
                    best_idx = jj
                end
            end
        end
        used[best_idx] = true
        out[:, ii] .= view(left_vectors, :, best_idx)
    end

    return out
end

"""
    _min_spectral_gaps(vals)

Returns minimum distance to any other eigenvalue for each index.
"""
function _min_spectral_gaps(vals::AbstractVector{Complex{BigFloat}})
    n = length(vals)
    gaps = fill(BigFloat(Inf), n)
    @inbounds for ii in 1:n
        g = BigFloat(Inf)
        for jj in 1:n
            if ii != jj
                d = abs(vals[ii] - vals[jj])
                if d < g
                    g = d
                end
            end
        end
        gaps[ii] = g
    end
    return gaps
end

"""
    _two_sided_refine_pair(A_big, A_adj_big, λ0, x0, y0, max_iter, tol1; damping=0.8, trust_radius=0.5, biorth_tol=1e-30)

Two-sided nonsymmetric refinement with right/left residual corrections.
"""
function _two_sided_refine_pair(
    A_big::AbstractMatrix{Complex{BigFloat}},
    A_adj_big::AbstractMatrix{Complex{BigFloat}},
    λ0::Complex{BigFloat},
    x0::AbstractVector{Complex{BigFloat}},
    y0::AbstractVector{Complex{BigFloat}},
    max_iter::Integer,
    tol1::BigFloat;
    damping::Float64=0.8,
    trust_radius::Float64=0.5,
    biorth_tol::Float64=1e-30,
)
    DT = BigFloat
    n = length(x0)

    x = copy(x0)
    y = copy(y0)
    λ = λ0

    damping_dt = DT(damping)
    trust_radius_dt = DT(trust_radius)
    biorth_tol_dt = DT(biorth_tol)

    mode = :iterative_maxiter
    warning = ""
    rr = BigFloat(Inf)
    rl = BigFloat(Inf)
    berr = BigFloat(Inf)
    bierr = BigFloat(Inf)
    cond_proxy = BigFloat(Inf)

    for iter in 1:max_iter
        nx = max(norm(x, Inf), eps(DT))
        x ./= nx

        β = dot(y, x)
        if abs(β) <= biorth_tol_dt
            mode = :biorthogonality_failed
            warning = "biorthogonality denominator collapsed"
            return λ, x, y, iter, false, mode, rr, rl, berr, bierr, cond_proxy, warning
        end
        y ./= β

        Ax = A_big * x
        Aty = A_adj_big * y

        rR = Ax .- λ .* x
        rL = Aty .- conj(λ) .* y

        rr = norm(rR, Inf) / max(abs(λ) * norm(x, Inf), eps(DT))
        rl = norm(rL, Inf) / max(abs(λ) * norm(y, Inf), eps(DT))

        berr = max(rr, rl)
        bierr = abs(dot(y, x) - one(Complex{BigFloat}))
        cond_proxy = norm(x, Inf) * norm(y, Inf)

        if berr <= tol1 && bierr <= sqrt(tol1)
            mode = :iterative
            return λ, x, y, iter, true, mode, rr, rl, berr, bierr, cond_proxy, warning
        end

        M = copy(A_big)
        M[diagind(M)] .-= λ
        s = argmax(abs.(x))
        M[:, s] .= -x

        N = copy(A_adj_big)
        N[diagind(N)] .-= conj(λ)
        t = argmax(abs.(y))
        N[:, t] .= -y

        δx = nothing
        δy = nothing
        try
            Fx = lu(M)
            Fy = lu(N)
            δx = -(Fx \ rR)
            δy = -(Fy \ rL)
        catch
            mode = :factorization_failed
            warning = "correction system factorization failed"
            return λ, x, y, iter, false, mode, rr, rl, berr, bierr, cond_proxy, warning
        end

        sx = norm(δx, Inf) / max(norm(x, Inf), eps(DT))
        sy = norm(δy, Inf) / max(norm(y, Inf), eps(DT))
        step = max(sx, sy)

        α = damping_dt
        if step > trust_radius_dt
            α *= trust_radius_dt / max(step, eps(DT))
        end

        x .+= α .* δx
        y .+= α .* δy

        Ax2 = A_big * x
        β2 = dot(y, x)
        if abs(β2) <= biorth_tol_dt
            mode = :biorthogonality_failed
            warning = "biorthogonality denominator collapsed during update"
            return λ, x, y, iter, false, mode, rr, rl, berr, bierr, cond_proxy, warning
        end
        λ = dot(y, Ax2) / β2
    end

    return λ, x, y, max_iter, false, mode, rr, rl, berr, bierr, cond_proxy, warning
end

"""
    EigenRefNonSym(A; Low_pres_type=Float64, Refinement_precision=256, Max_iter=40, tol_fact=1, cluster_gap_factor=1e3, damping=0.8, trust_radius=0.5, return_status=false, show_progress=false)

Refines eigenpairs of a general real nonsymmetric matrix using two-sided iterative refinement.
Returns `(vals, right_vecs, left_vecs)` and optionally `status`.

If `return_status=true`, `status` includes:
- `mode_used`, `iterations_used`, `converged`
- `right_residual`, `left_residual`, `backward_error`
- `biorthogonality_error`, `condition_proxy`, `clustered`, `warnings`
- `restrictions` (scope of guaranteed convergence)
"""
function EigenRefNonSym(
    A::Array{T,2};
    Low_pres_type::Type=Float64,
    Refinement_precision::Integer=256,
    Max_iter::Integer=40,
    tol_fact::Integer=1,
    cluster_gap_factor::Float64=1e3,
    damping::Float64=0.8,
    trust_radius::Float64=0.5,
    return_status::Bool=false,
    show_progress::Bool=false,
) where {T <: Union{Float32, Float64, BigFloat}}
    A_low = Low_pres_type.(A)
    er = eigen(A_low)
    el = eigen(adjoint(A_low))

    vals0 = ComplexF64.(er.values)
    vecs0 = ComplexF64.(er.vectors)
    left0 = _match_left_vectors(vals0, ComplexF64.(el.values), ComplexF64.(el.vectors))

    return EigenRefNonSym(
        A,
        vals0,
        vecs0;
        left_vecs0=left0,
        Refinement_precision=Refinement_precision,
        Max_iter=Max_iter,
        tol_fact=tol_fact,
        cluster_gap_factor=cluster_gap_factor,
        damping=damping,
        trust_radius=trust_radius,
        return_status=return_status,
        show_progress=show_progress,
    )
end

"""
    EigenRefNonSym(A, vals0, vecs0; left_vecs0=nothing, ...)

Refines from user-provided initial eigenvalues/eigenvectors.

When `return_status=true`, status diagnostics mirror the fields documented in
the primary `EigenRefNonSym(A; ...)` entry point.
"""
function EigenRefNonSym(
    A::Array{T,2},
    vals0::AbstractVector,
    vecs0::AbstractMatrix;
    left_vecs0=nothing,
    Refinement_precision::Integer=256,
    Max_iter::Integer=40,
    tol_fact::Integer=1,
    cluster_gap_factor::Float64=1e3,
    damping::Float64=0.8,
    trust_radius::Float64=0.5,
    return_status::Bool=false,
    show_progress::Bool=false,
) where {T <: Union{Float32, Float64, BigFloat}}
    n = size(A, 1)
    if size(A, 1) != size(A, 2)
        throw(ArgumentError("A must be square"))
    end
    if size(vecs0, 1) != n || size(vecs0, 2) != length(vals0)
        throw(ArgumentError("Initial vecs0 dimensions must match A and vals0"))
    end

    return setprecision(BigFloat, Refinement_precision) do
        DT = BigFloat
        tol1 = tol_fact * eps(DT)

        A_big = Complex{DT}.(A)
        A_adj_big = adjoint(A_big)

        vals_init = Complex{DT}.(vals0)
        right_init = Complex{DT}.(vecs0)

        left_init = if isnothing(left_vecs0)
            el = eigen(adjoint(ComplexF64.(A)))
            Complex{DT}.(_match_left_vectors(ComplexF64.(vals0), ComplexF64.(el.values), ComplexF64.(el.vectors)))
        else
            if size(left_vecs0, 1) != n || size(left_vecs0, 2) != length(vals0)
                throw(ArgumentError("left_vecs0 dimensions must match A and vals0"))
            end
            Complex{DT}.(left_vecs0)
        end

        m = length(vals_init)
        vals_out = Vector{Complex{DT}}(undef, m)
        right_out = Matrix{Complex{DT}}(undef, n, m)
        left_out = Matrix{Complex{DT}}(undef, n, m)

        mode_used = return_status ? Vector{Symbol}(undef, m) : Symbol[]
        iterations_used = return_status ? zeros(Int64, m) : Int64[]
        converged = return_status ? falses(m) : BitVector()
        right_residual = return_status ? fill(BigFloat(Inf), m) : BigFloat[]
        left_residual = return_status ? fill(BigFloat(Inf), m) : BigFloat[]
        backward_error = return_status ? fill(BigFloat(Inf), m) : BigFloat[]
        biorth_error = return_status ? fill(BigFloat(Inf), m) : BigFloat[]
        condition_proxy = return_status ? fill(BigFloat(Inf), m) : BigFloat[]
        warnings = return_status ? fill("", m) : String[]
        progress_state = _progress_init(m; show=show_progress, label="Nonsymmetric refinement")

        gaps = _min_spectral_gaps(vals_init)
        cluster_tol = DT(cluster_gap_factor) * sqrt(eps(DT))
        clustered = gaps .<= cluster_tol .* max.(abs.(vals_init), one(DT))

        # Schur fallback used for clustered/near-defective pairs.
        S = schur(ComplexF64.(A))
        schur_vals = Complex{DT}.(S.values)
        schur_vecs = Complex{DT}.(S.Z)

        SL = schur(adjoint(ComplexF64.(A)))
        schur_left = Complex{DT}.(_match_left_vectors(ComplexF64.(S.values), ComplexF64.(SL.values), ComplexF64.(SL.vectors)))

        function _refine_nonsym_index!(ii)
            λ0 = vals_init[ii]
            x0 = copy(view(right_init, :, ii))
            y0 = copy(view(left_init, :, ii))

            if clustered[ii]
                idx = argmin(abs.(schur_vals .- λ0))
                x = copy(view(schur_vecs, :, idx))
                y = copy(view(schur_left, :, idx))
                β = dot(y, x)
                if abs(β) > eps(DT)
                    y ./= β
                end
                λ = schur_vals[idx]

                vals_out[ii] = λ
                right_out[:, ii] .= x
                left_out[:, ii] .= y

                if return_status
                    mode_used[ii] = :schur_block_fallback
                    iterations_used[ii] = 0
                    converged[ii] = false
                    rR = A_big * x .- λ .* x
                    rL = A_adj_big * y .- conj(λ) .* y
                    rr = norm(rR, Inf) / max(abs(λ) * norm(x, Inf), eps(DT))
                    rl = norm(rL, Inf) / max(abs(λ) * norm(y, Inf), eps(DT))
                    right_residual[ii] = rr
                    left_residual[ii] = rl
                    backward_error[ii] = max(
                        rr,
                        rl,
                    )
                    biorth_error[ii] = abs(dot(y, x) - one(Complex{DT}))
                    condition_proxy[ii] = norm(x, Inf) * norm(y, Inf)
                    warnings[ii] = "clustered spectrum; using Schur fallback"
                end
                _progress_update!(progress_state)
                return nothing
            end

            λ, x, y, iter, ok, mode, rr, rl, berr, bierr, condp, warn = _two_sided_refine_pair(
                A_big,
                A_adj_big,
                λ0,
                x0,
                y0,
                Max_iter,
                tol1;
                damping=damping,
                trust_radius=trust_radius,
            )

            vals_out[ii] = λ
            right_out[:, ii] .= x
            left_out[:, ii] .= y

            if return_status
                mode_used[ii] = mode
                iterations_used[ii] = iter
                converged[ii] = ok
                right_residual[ii] = rr
                left_residual[ii] = rl
                backward_error[ii] = berr
                biorth_error[ii] = bierr
                condition_proxy[ii] = condp
                warnings[ii] = warn
            end

            _progress_update!(progress_state)
            return nothing
        end

        let left_init = left_init, schur_left = schur_left
            OhMyThreads.@tasks for ii in 1:m
                OhMyThreads.@set scheduler = :dynamic
                OhMyThreads.@set ntasks = Threads.nthreads()
                _refine_nonsym_index!(ii)
            end
        end

        _progress_finish!(progress_state)

        if return_status
            status = (
                mode_used = mode_used,
                iterations_used = iterations_used,
                converged = converged,
                right_residual = right_residual,
                left_residual = left_residual,
                backward_error = backward_error,
                biorthogonality_error = biorth_error,
                condition_proxy = condition_proxy,
                clustered = clustered,
                warnings = warnings,
                restrictions = "Guaranteed convergence requires a diagonalizable target pair and non-clustered spectrum; strongly non-normal cases are best-effort.",
            )
            return vals_out, right_out, left_out, status
        end

        return vals_out, right_out, left_out
    end
end
