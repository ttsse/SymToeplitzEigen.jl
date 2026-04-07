"""
    _is_symmetric_toeplitz(A)

Checks whether matrix `A` is approximately symmetric Toeplitz.
"""
function _is_symmetric_toeplitz(A :: AbstractMatrix{T}) where {T <: Union{Float32, Float64}}
    n1, n2 = size(A)
    if n1 != n2
        return false
    end

    atol = 100 * eps(T)
    @inbounds begin
        for jj in 1:n1
            for kk in 1:n1
                ref = A[abs(kk - jj) + 1, 1]
                if !isapprox(A[kk, jj], ref; atol=atol, rtol=0)
                    return false
                end
            end
        end
    end
    return true
end

"""
    _refine_pair_robust(A, λ0, x0, precision_bits, max_iter, tol_fact; toeplitz_col = nothing)

Refines one eigenpair in robust (BigFloat solve) mode and returns
`(λ, x, iter, converged)`.
"""
function _refine_pair_robust(
    A :: Array{T, 2},
    λ0 :: BigFloat,
    x0 :: AbstractVector{BigFloat},
    precision_bits :: Integer,
    max_iter :: Integer,
    tol_fact :: Integer,
    ;
    toeplitz_col = nothing,
) where T
    return setprecision(BigFloat, precision_bits) do
        DT = BigFloat
        n = length(x0)
        use_toeplitz = !isnothing(toeplitz_col)

        A_big = use_toeplitz ? nothing : DT.(A)
        c_big = use_toeplitz ? DT.(toeplitz_col) : DT[]
        Bbig = DT.(zeros(Int64, n, n))
        r = DT.(zeros(Int64, n))
        y = DT.(zeros(Int64, n))
        yp = DT.(zeros(Int64, n))
        x = DT.(zeros(Int64, n))

        my_vec_setvalue!(n, x, x0)
        λ = DT(λ0)
        s = normalize_Infnorm!(x)

        tol1 = tol_fact * eps(DT)
        err_norm = DT(Inf)
        iter = 0

        # Initial factorization corresponds to the first iteration state.
        if use_toeplitz
            my_fill_sym_toeplitz!(n, Bbig, c_big)
        else
            my_mat_setvalue!(n, Bbig, A_big)
        end
        my_add_diag_elements!(n, Bbig, -λ)
        @inbounds begin
            for kk in 1:n
                setvalue_BF!(Bbig[kk, s], x[kk])
                neg!(Bbig[kk, s])
            end
        end
        FB = lu!(Bbig)

        while err_norm > tol1 && iter < max_iter
            iter += 1

            if use_toeplitz
                my_neg_sym_toeplitz_mat_vec!(n, r, c_big, x)
            else
                my_neg_mat_vec_mul!(n, r, A_big, x)
            end
            my_add_scaled_vec!(n, r, λ, x)

            if iter == 2
                if use_toeplitz
                    my_fill_sym_toeplitz!(n, Bbig, c_big)
                else
                    my_mat_setvalue!(n, Bbig, A_big)
                end
                my_add_diag_elements!(n, Bbig, -λ)
                @inbounds begin
                    for kk in 1:n
                        setvalue_BF!(Bbig[kk, s], x[kk])
                        neg!(Bbig[kk, s])
                    end
                end
                FB = lu!(Bbig)
            end

            ldiv!(y, (FB::LinearAlgebra.LU), r)
            my_vec_setvalue!(n, yp, y)

            ys = y[s]
            yp[s] = zero(DT)

            my_add_vec!(n, x, yp)
            λ += ys

            r_norm = my_maxabs(n, r)
            x_norm = my_maxabs(n, x)
            err_norm = r_norm / max(abs(λ) * x_norm, eps(DT))
        end

        return λ, x, iter, err_norm <= tol1
    end
end

"""
    Refinement(Tn, vals, vecs, Refinement_precision, Max_iter, tol_fact)

Improves the accuracy of the eigenvalues `vals` (and vectors `vecs`) of `Tn` up to a precision given by `tol_fact` * eps(BigFloat) where BigFloat will use `Refinement_precision` bits, with a maximum of `Max_iter` iterations.
"""
function Refinement(
    A :: Array{T, 2},
    Vals :: Array{T, 1},
    Vecs :: Array{T, 2},
    Refinement_precision :: Integer,
    Max_iter :: Integer,
    tol_fact :: Integer;
    solve_mode :: Symbol = :fast,
    stall_ratio :: Float64 = 0.9,
    stall_iters :: Integer = 3,
    return_status :: Bool = false,
    adaptive_precision_escalation :: Bool = false,
    escalation_precision :: Integer = 0,
    escalation_extra_iter :: Integer = 10,
    toeplitz_kernel :: Symbol = :auto,
    toeplitz_column = nothing,
    toeplitz_auto_threshold :: Integer = 96,
    show_progress :: Bool = false,
) where T
    if !(solve_mode in (:fast, :robust, :adaptive))
        throw(ArgumentError("Invalid solve_mode=$(solve_mode). Expected one of :fast, :robust, or :adaptive."))
    end
    if stall_iters < 1
        throw(ArgumentError("stall_iters must be >= 1."))
    end
    if !(0.0 < stall_ratio <= 1.0)
        throw(ArgumentError("stall_ratio must satisfy 0 < stall_ratio <= 1."))
    end
    if !(toeplitz_kernel in (:dense, :structured, :auto, :fft))
        throw(ArgumentError("Invalid toeplitz_kernel=$(toeplitz_kernel). Expected :dense, :structured, :auto, or :fft."))
    end
    if escalation_extra_iter < 1
        throw(ArgumentError("escalation_extra_iter must be >= 1."))
    end
    if toeplitz_auto_threshold < 1
        throw(ArgumentError("toeplitz_auto_threshold must be >= 1."))
    end

    target_escalation_precision = escalation_precision > 0 ? escalation_precision : 2 * Refinement_precision
    if adaptive_precision_escalation && target_escalation_precision <= Refinement_precision
        throw(ArgumentError("escalation_precision must be larger than Refinement_precision when adaptive_precision_escalation=true."))
    end

    setprecision(BigFloat, Refinement_precision)
    DT = BigFloat
    stall_ratio_dt = DT(stall_ratio)

    m = length(Vals)
    n = size(Vecs)[1]
    refinedVals = DT.(zeros(Int64, m))
    refinedVecs = DT.(zeros(Int64, n, m))

    effective_toeplitz_kernel = :dense
    c_src = nothing
    if toeplitz_kernel == :structured
        if isnothing(toeplitz_column)
            if !_is_symmetric_toeplitz(A)
                throw(ArgumentError("toeplitz_kernel=:structured requires a symmetric Toeplitz matrix input, or provide toeplitz_column explicitly."))
            end
            c_src = @view A[:, 1]
        else
            c_src = toeplitz_column
        end
        effective_toeplitz_kernel = :structured
    elseif toeplitz_kernel == :fft
        if isnothing(toeplitz_column)
            if !_is_symmetric_toeplitz(A)
                throw(ArgumentError("toeplitz_kernel=:fft requires a symmetric Toeplitz matrix input, or provide toeplitz_column explicitly."))
            end
            c_src = @view A[:, 1]
        else
            c_src = toeplitz_column
        end
        effective_toeplitz_kernel = :fft
    elseif toeplitz_kernel == :auto
        if n >= toeplitz_auto_threshold
            if isnothing(toeplitz_column)
                if _is_symmetric_toeplitz(A)
                    c_src = @view A[:, 1]
                    effective_toeplitz_kernel = :structured
                end
            else
                c_src = toeplitz_column
                effective_toeplitz_kernel = :structured
            end
        end
    end

    resolved_toeplitz_kernel =
        (effective_toeplitz_kernel == :fft && solve_mode != :fast) ? :structured : effective_toeplitz_kernel

    use_toeplitz = resolved_toeplitz_kernel in (:structured, :fft)

    A_big = use_toeplitz ? nothing : DT.(A)
    A_f64 = use_toeplitz ? nothing : Float64.(A)
    c_big = use_toeplitz ? DT.(c_src) : DT[]
    c_f64 = use_toeplitz ? Float64.(c_src) : Float64[]

    tol1 = tol_fact * eps(DT)

    mode_used = return_status ? Vector{Symbol}(undef, m) : Symbol[]
    iterations_used = return_status ? zeros(Int64, m) : Int64[]
    converged = return_status ? falses(m) : BitVector()
    promoted_iteration = return_status ? zeros(Int64, m) : Int64[]
    precision_bits_used = return_status ? fill(Int64(Refinement_precision), m) : Int64[]
    precision_escalated = return_status ? falses(m) : BitVector()
    progress_state = _progress_init(m; show=show_progress, label="Symmetric refinement")

    function _refine_pair_index!(
        ii,
        Btmp,
        r,
        y,
        yp,
        r64,
        y64,
        x64,
        Bbig_ref,
        fft_cache,
    )
        λ = DT(Vals[ii])
        x = @view refinedVecs[:, ii]
        my_vec_setvalue_prom!(n, x, @view Vecs[:, ii])

        s = normalize_Infnorm!(x)

        iter = 0
        stalled_count = 0
        mode_promoted = false
        promoted_iter = 0

        err_norm = DT(Inf)
        prev_err_norm = DT(Inf)
        current_mode = solve_mode == :adaptive ? :fast : solve_mode
        FB_fast = nothing
        FB_robust = nothing
        rebuild_fast = true
        rebuild_robust = true

        while err_norm > tol1 && iter < Max_iter
            iter += 1

            # Compute residual r = (λI - A)x in BigFloat without touching shared matrix state.
            if use_toeplitz
                if resolved_toeplitz_kernel == :fft
                    my_neg_sym_toeplitz_mat_vec_fft!(n, r, c_f64, x, x64, fft_cache)
                else
                    my_neg_sym_toeplitz_mat_vec!(n, r, c_big, x)
                end
            else
                my_neg_mat_vec_mul!(n, r, A_big, x)
            end
            my_add_scaled_vec!(n, r, λ, x)

            if current_mode == :fast
                # Since Btmp is lower precision, after the second iteration no need to refresh each loop.
                if iter < 3 || rebuild_fast
                    if use_toeplitz
                        my_fill_sym_toeplitz!(n, Btmp, c_f64)
                    else
                        copyto!(Btmp, A_f64)
                    end
                    Btmp[diagind(Btmp)] .-= Float64(λ)
                    my_vec_demote!(n, x64, x)
                    @inbounds begin
                        for kk in 1:n
                            Btmp[kk, s] = -x64[kk]
                        end
                    end

                    FB_fast = lu!(Btmp)
                    rebuild_fast = false
                end

                # Solve low-precision system then promote correction to BigFloat.
                my_vec_demote!(n, r64, r)
                ldiv!(y64, (FB_fast::LinearAlgebra.LU), r64)
                my_vec_setvalue_prom!(n, y, y64)
            else
                Bbig = Bbig_ref[]
                if isnothing(Bbig)
                    Bbig = DT.(zeros(Int64, n, n))
                    Bbig_ref[] = Bbig
                end
                Bbig_mat = Bbig::Matrix{DT}

                # Robust mode factors and solves directly in BigFloat.
                if iter < 3 || rebuild_robust
                    if use_toeplitz
                        my_fill_sym_toeplitz!(n, Bbig_mat, c_big)
                    else
                        my_mat_setvalue!(n, Bbig_mat, A_big)
                    end
                    my_add_diag_elements!(n, Bbig_mat, -λ)
                    @inbounds begin
                        for kk in 1:n
                            setvalue_BF!(Bbig_mat[kk, s], x[kk])
                            neg!(Bbig_mat[kk, s])
                        end
                    end

                    FB_robust = lu!(Bbig_mat)
                    rebuild_robust = false
                end

                ldiv!(y, (FB_robust::LinearAlgebra.LU), r)
            end

            my_vec_setvalue!(n, yp, y)

            ys = y[s]
            yp[s] = zero(DT)

            # Update
            my_add_vec!(n, x, yp)
            λ += ys

            # Deterministic relative residual estimate.
            r_norm = my_maxabs(n, r)
            x_norm = my_maxabs(n, x)
            err_norm = r_norm / max(abs(λ) * x_norm, eps(DT))

            if solve_mode == :adaptive && current_mode == :fast && iter > 1
                reduction = err_norm / max(prev_err_norm, eps(DT))
                if reduction >= stall_ratio_dt
                    stalled_count += 1
                else
                    stalled_count = 0
                end

                if stalled_count >= stall_iters
                    current_mode = :robust
                    mode_promoted = true
                    promoted_iter = iter
                    stalled_count = 0
                    rebuild_robust = true
                end
            end

            prev_err_norm = err_norm
        end

        converged_pair = err_norm <= tol1
        did_precision_escalate = false
        bits_used = Refinement_precision

        if solve_mode == :adaptive && adaptive_precision_escalation && !converged_pair
            λ_hi, x_hi, iter_hi, converged_hi = _refine_pair_robust(
                A,
                λ,
                x,
                target_escalation_precision,
                escalation_extra_iter,
                tol_fact,
                toeplitz_col = use_toeplitz ? c_big : nothing,
            )

            λ = λ_hi
            @views refinedVecs[:, ii] .= x_hi
            iter += iter_hi
            converged_pair = converged_hi
            did_precision_escalate = true
            bits_used = target_escalation_precision
        end

        refinedVals[ii] = λ

        if return_status
            iterations_used[ii] = iter
            converged[ii] = converged_pair
            promoted_iteration[ii] = promoted_iter
            precision_bits_used[ii] = bits_used
            precision_escalated[ii] = did_precision_escalate

            if solve_mode == :adaptive
                if did_precision_escalate
                    mode_used[ii] = converged_pair ? :adaptive_precision_escalated : :adaptive_precision_escalated_maxiter
                elseif mode_promoted
                    mode_used[ii] = converged_pair ? :adaptive_promoted : :adaptive_promoted_maxiter
                else
                    mode_used[ii] = converged_pair ? :adaptive_fast : :adaptive_fast_maxiter
                end
            elseif solve_mode == :fast
                mode_used[ii] = converged_pair ? :fast : :fast_maxiter
            else
                mode_used[ii] = converged_pair ? :robust : :robust_maxiter
            end
        end

        _progress_update!(progress_state)
        return nothing
    end

    OhMyThreads.@tasks for ii in 1:m
        @set scheduler = :dynamic
        @set ntasks = Threads.nthreads()
        @local ws = (
            Btmp = Matrix{Float64}(undef, n, n),
            r = DT.(zeros(Int64, n)),
            y = DT.(zeros(Int64, n)),
            yp = DT.(zeros(Int64, n)),
            r64 = Vector{Float64}(undef, n),
            y64 = Vector{Float64}(undef, n),
            x64 = Vector{Float64}(undef, n),
            Bbig_ref = Ref{Union{Nothing, Matrix{DT}}}(nothing),
            fft_cache = Ref{Any}(nothing),
        )

        _refine_pair_index!(ii, ws.Btmp, ws.r, ws.y, ws.yp, ws.r64, ws.y64, ws.x64, ws.Bbig_ref, ws.fft_cache)
    end

    _progress_finish!(progress_state)

    if return_status
        status = (
            mode_used = mode_used,
            iterations_used = iterations_used,
            converged = converged,
            promoted_iteration = promoted_iteration,
            precision_bits_used = precision_bits_used,
            precision_escalated = precision_escalated,
            toeplitz_kernel_used = resolved_toeplitz_kernel,
        )
        return refinedVals, refinedVecs, status
    end

    return refinedVals, refinedVecs
end
