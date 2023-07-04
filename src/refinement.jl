"""
    Refinement(Tn, vals, vecs, Refinement_precision, Max_iter, tol_fact)

Improves the accuracy of the eigenvalues `vals` (and vectors `vecs`) of `Tn` up to a precision given by `tol_fact` * eps(BigFloat) where BigFloat will use `Refinement_precision` bits, with a maximum of `Max_iter` iterations.
"""
function Refinement(A :: Array{T, 2}, Vals :: Array{T, 1}, Vecs :: Array{T, 2}, Refinement_precision :: Integer, Max_iter :: Integer, tol_fact :: Integer) where T
    setprecision(BigFloat, Refinement_precision)
    DT = BigFloat

    m = length(Vals)
    n = size(Vecs)[1]
    refinedVals = DT.(zeros(Int64,n))
    refinedVecs = DT.(zeros(Int64,n,n))

    tol1 = tol_fact * eps(DT)

    # Prepare for multi-threading
    num_threads = Threads.nthreads()
    len, rem = divrem(m, num_threads)
    
    Threads.@threads for pp = 1:num_threads

        f = 1 + ((pp-1) * len)
        l = f + len - 1

        if rem > 0
            if pp <= rem
                f = f + (pp-1)
                l = l + pp
            else
                f = f + rem
                l = l + rem
            end
        end

        

        # Pre-allocate
        B = DT.(A)
        Btmp = zeros(T, n, n)
        r = DT.(zeros(Int64,n))
        y = DT.(zeros(Int64,n))
        ys = DT(0)
        yp = DT.(zeros(Int64,n))
        

        for ii in f:l

            λ = DT(Vals[ii])
            x = DT.(@view Vecs[:, ii])

            s = normalize_Infnorm!(x)

            iter = 0

            err_norm = DT(Inf)
            FB = nothing

            while err_norm > tol1 && iter < Max_iter
                iter += 1
                
                # We compute `r = (B - λI)x` so first we will remove λ from the diagonal and then re-add it after 
                # the computation so we can reuse the same B in all iterations without needing to allocate
                my_add_diag_elements!(n, B, -1*λ)
                
                # Compute Residual
                my_neg_mat_vec_mul!(n, r, B, x)

                my_add_diag_elements!(n, B, λ)

                # Since Btmp is of a less accurate data type, after the second iteration no need to do this update
                if iter < 3 
                    Btmp .= A
                    Btmp[diagind(Btmp)] .-= T.(λ) 
                    Btmp[:,s] .= x
                    Btmp[:,s] *= -1
                    FB = lu!(Btmp)
                end

                # Solve Btmp*y = r
                my_vec_setvalue_prom!(n, y, FB \ T.(r))
                
                my_vec_setvalue!(n, yp, y)

                ys = y[s]
                yp[s] = zero(DT)
                
                # Update
                x += yp
                λ += ys

                # Cheap estimate error
                ee = rand(1:n)
                err_norm = abs((dot(view(A, ee, :), x) - λ*x[ee])/(λ*x[ee]))
            end
            refinedVals[ii] = λ
            refinedVecs[:,ii] .= x

        end
    end
    return refinedVals, refinedVecs
end