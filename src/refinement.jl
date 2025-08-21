"""
    Refinement(Tn, vals, vecs, Refinement_precision, Max_iter, tol_fact; track_errors=false)

Improves the accuracy of the eigenvalues `vals` (and vectors `vecs`) of `Tn` up to a precision
given by `tol_fact * eps(BigFloat)` where BigFloat will use `Refinement_precision` bits,
with a maximum of `Max_iter` iterations.

If `track_errors=true`, also returns a dictionary of error histories for each eigenvalue.
"""
function Refinement(A::Array{T,2}, Vals::Array{T,1}, Vecs::Array{T,2},
                    Refinement_precision::Integer, Max_iter::Integer, tol_fact::Integer;
                    track_errors::Bool=false) where T
    setprecision(BigFloat, Refinement_precision)
    DT = BigFloat

    m = length(Vals)
    n = size(Vecs, 1)
    refinedVals = DT.(zeros(Int64,n))
    refinedVecs = DT.(zeros(Int64,n,n))

    tol1 = tol_fact * eps(DT)

    error_history = track_errors ? Dict{Int, Vector{Float64}}() : nothing

    Threads.@threads for ii in 1:m
        λ = DT(Vals[ii])
        x = DT.(@view Vecs[:, ii])
        s = normalize_Infnorm!(x)

        err_norm = DT(Inf)
        iter = 0
        FB = nothing

        local_errors = track_errors ? Float64[] : nothing

        B = PermutedDimsArray(DT.(A), (2,1))
        Btmp = zeros(Float64, n, n)
        r = DT.(zeros(Int64,n))
        y = DT.(zeros(Int64,n))
        ys = DT(0)
        yp = DT.(zeros(Int64,n))

        while err_norm > tol1 && iter < Max_iter
            iter += 1

            my_add_diag_elements!(n, B, -1*λ)
            my_neg_mat_vec_mul!(n, r, B, x)
            my_add_diag_elements!(n, B, λ)

            if iter < 3
                Btmp .= A
                Btmp[diagind(Btmp)] .-= Float64(λ)
                Btmp[:,s] .= x
                Btmp[:,s] *= -1
                FB = lu!(Btmp)
            end

            my_vec_setvalue_prom!(n, y, FB \ Float64.(r))
            my_vec_setvalue!(n, yp, y)

            ys = y[s]
            yp[s] = zero(DT)

            x += yp
            λ += ys

            ee = rand(1:n)
            err_norm = abs((dot(view(A, ee, :), x) - λ*x[ee])/(λ*x[ee]))

            if track_errors
                push!(local_errors, Float64(err_norm))
            end
        end

        refinedVals[ii] = λ
        refinedVecs[:,ii] .= x
        if track_errors
            error_history[ii] = local_errors
        end
    end

    return track_errors ? (refinedVals, refinedVecs, error_history) : (refinedVals, refinedVecs)
end
