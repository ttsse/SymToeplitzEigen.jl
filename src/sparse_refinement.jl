function _sparse_bandwidth_ratio(A::SparseMatrixCSC{T, Ti}) where {T, Ti <: Integer}
    n = size(A, 1)
    n == 0 && return 0.0

    max_band = 0
    rows = rowvals(A)
    @inbounds begin
        for jj in 1:n
            for ptr in nzrange(A, jj)
                ii = rows[ptr]
                d = abs(ii - jj)
                if d > max_band
                    max_band = d
                end
            end
        end
    end

    return (2 * max_band + 1) / n
end

function RefinementSparseFast(
    A::SparseMatrixCSC{Float64, Ti},
    Vals::AbstractVector,
    Vecs::AbstractMatrix,
    Refinement_precision::Integer,
    Max_iter::Integer,
    tol_fact::Integer;
    return_status::Bool=false,
    show_progress::Bool=false,
) where {Ti <: Integer}
    setprecision(BigFloat, Refinement_precision)
    DT = BigFloat

    m = length(Vals)
    n = size(A, 1)
    refinedVals = DT.(zeros(Int64, m))
    refinedVecs = DT.(zeros(Int64, n, m))

    tol1 = tol_fact * eps(DT)

    mode_used = return_status ? Vector{Symbol}(undef, m) : Symbol[]
    iterations_used = return_status ? zeros(Int64, m) : Int64[]
    converged = return_status ? falses(m) : BitVector()
    promoted_iteration = return_status ? zeros(Int64, m) : Int64[]
    precision_bits_used = return_status ? fill(Int64(Refinement_precision), m) : Int64[]
    precision_escalated = return_status ? falses(m) : BitVector()

    progress_state = _progress_init(m; show=show_progress, label="Sparse refinement")

    function _refine_sparse_pair!(ii, r, y, yp, r64, y64, x64)
        λ = DT(Vals[ii])
        x = @view refinedVecs[:, ii]
        my_vec_setvalue_prom!(n, x, @view Vecs[:, ii])

        s = normalize_Infnorm!(x)
        iter = 0
        err_norm = DT(Inf)
        FB_fast = nothing
        rebuild_fast = true

        while err_norm > tol1 && iter < Max_iter
            iter += 1

            my_neg_sparse_mat_vec_mul!(n, r, A, x, x64, r64)
            my_add_scaled_vec!(n, r, λ, x)

            if iter < 3 || rebuild_fast
                my_vec_demote!(n, x64, x)
                Bsp = copy(A)
                λ64 = Float64(λ)

                @inbounds begin
                    for kk in 1:n
                        Bsp[kk, kk] -= λ64
                        Bsp[kk, s] = -x64[kk]
                    end
                end

                FB_fast = lu(Bsp)
                rebuild_fast = false
            end

            my_vec_demote!(n, r64, r)
            ldiv!(y64, FB_fast, r64)
            my_vec_setvalue_prom!(n, y, y64)
            my_vec_setvalue!(n, yp, y)

            ys = y[s]
            yp[s] = zero(DT)

            my_add_vec!(n, x, yp)
            λ += ys

            r_norm = my_maxabs(n, r)
            x_norm = my_maxabs(n, x)
            err_norm = r_norm / max(abs(λ) * x_norm, eps(DT))
        end

        refinedVals[ii] = λ

        if return_status
            iterations_used[ii] = iter
            converged[ii] = err_norm <= tol1
            promoted_iteration[ii] = 0
            precision_bits_used[ii] = Refinement_precision
            precision_escalated[ii] = false
            mode_used[ii] = converged[ii] ? :fast : :fast_maxiter
        end

        _progress_update!(progress_state)
        return nothing
    end

    OhMyThreads.@tasks for ii in 1:m
        @set scheduler = :dynamic
        @set ntasks = Threads.nthreads()
        @local ws = (
            r = DT.(zeros(Int64, n)),
            y = DT.(zeros(Int64, n)),
            yp = DT.(zeros(Int64, n)),
            r64 = Vector{Float64}(undef, n),
            y64 = Vector{Float64}(undef, n),
            x64 = Vector{Float64}(undef, n),
        )

        _refine_sparse_pair!(ii, ws.r, ws.y, ws.yp, ws.r64, ws.y64, ws.x64)
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
            toeplitz_kernel_used = :dense,
            sparse_solver_used = :sparse,
        )
        return refinedVals, refinedVecs, status
    end

    return refinedVals, refinedVecs
end

"""
    EigenRef(A::SparseMatrixCSC; ...)

Sparse-matrix overload for symmetric inputs. By default, this performs the same
full-spectrum workflow as dense matrices and auto-selects a sparse correction
path when solve_mode=:fast and sparsity is favorable.

Sparse inputs always initialize and refine the full spectrum.
"""
function EigenRef(
    A::SparseMatrixCSC{T, Ti};
    Low_pres_type::Type=Float64,
    Refinement_precision::Integer=256,
    Max_iter::Integer=100,
    tol_fact::Integer=1,
    solve_mode::Symbol=:fast,
    stall_ratio::Float64=0.9,
    stall_iters::Integer=3,
    return_status::Bool=false,
    adaptive_precision_escalation::Bool=false,
    escalation_precision::Integer=0,
    escalation_extra_iter::Integer=10,
    show_progress::Bool=false,
    sparse_solver::Symbol=:auto,
    sparse_density_threshold::Float64=0.35,
    sparse_bandwidth_threshold::Float64=0.5,
    sparse_min_n::Integer=128,
) where {T <: Union{Float32, Float64}, Ti <: Integer}
    if !issymmetric(A)
        throw(ArgumentError("Given sparse matrix must be symmetric"))
    end
    if sparse_density_threshold <= 0 || sparse_density_threshold > 1
        throw(ArgumentError("sparse_density_threshold must satisfy 0 < threshold <= 1."))
    end
    if sparse_bandwidth_threshold <= 0 || sparse_bandwidth_threshold > 1
        throw(ArgumentError("sparse_bandwidth_threshold must satisfy 0 < threshold <= 1."))
    end
    if sparse_min_n < 1
        throw(ArgumentError("sparse_min_n must be >= 1."))
    end
    if !(sparse_solver in (:auto, :dense, :sparse))
        throw(ArgumentError("Invalid sparse_solver=$(sparse_solver). Expected :auto, :dense, or :sparse."))
    end

    n = size(A, 1)

    vals0, vecs0 = eigen(Symmetric(Matrix(Low_pres_type.(A))))

    vals = Float64.(vals0)
    vecs = Float64.(vecs0)

    density = nnz(A) / (n * n)
    bandwidth_ratio = _sparse_bandwidth_ratio(A)
    effective_sparse_solver = sparse_solver
    if sparse_solver == :auto
        use_sparse = solve_mode == :fast && n >= sparse_min_n && (
            density <= sparse_density_threshold || bandwidth_ratio <= sparse_bandwidth_threshold
        )
        effective_sparse_solver = use_sparse ? :sparse : :dense
    end

    if effective_sparse_solver == :sparse && solve_mode == :fast
        A64 = Float64.(A)
        if return_status
            nvals, nvecs, status = RefinementSparseFast(
                A64,
                vals,
                vecs,
                Refinement_precision,
                Max_iter,
                tol_fact;
                return_status=true,
                show_progress=show_progress,
            )

            status_with_meta = merge(
                status,
                (
                    sparse_solver_used=:sparse,
                    sparse_density=density,
                    sparse_bandwidth_ratio=bandwidth_ratio,
                ),
            )
            return nvals, nvecs, status_with_meta
        end

        return RefinementSparseFast(
            A64,
            vals,
            vecs,
            Refinement_precision,
            Max_iter,
            tol_fact;
            return_status=false,
            show_progress=show_progress,
        )
    end

    A_dense = Matrix{Float64}(A)
    if return_status
        nvals, nvecs, status = Refinement(
            A_dense,
            vals,
            vecs,
            Refinement_precision,
            Max_iter,
            tol_fact;
            solve_mode=solve_mode,
            stall_ratio=stall_ratio,
            stall_iters=stall_iters,
            return_status=true,
            adaptive_precision_escalation=adaptive_precision_escalation,
            escalation_precision=escalation_precision,
            escalation_extra_iter=escalation_extra_iter,
            toeplitz_kernel=:dense,
            show_progress=show_progress,
        )

        status_with_meta = merge(
            status,
            (
                sparse_solver_used=:dense,
                sparse_density=density,
                sparse_bandwidth_ratio=bandwidth_ratio,
            ),
        )

        return nvals, nvecs, status_with_meta
    end

    return Refinement(
        A_dense,
        vals,
        vecs,
        Refinement_precision,
        Max_iter,
        tol_fact;
        solve_mode=solve_mode,
        stall_ratio=stall_ratio,
        stall_iters=stall_iters,
        return_status=false,
        adaptive_precision_escalation=adaptive_precision_escalation,
        escalation_precision=escalation_precision,
        escalation_extra_iter=escalation_extra_iter,
        toeplitz_kernel=:dense,
        show_progress=show_progress,
    )
end

EigenRefSparse(A::SparseMatrixCSC; kwargs...) = EigenRef(A; kwargs...)
