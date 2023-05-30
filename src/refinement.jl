using IterativeRefinement
using LinearAlgebra

"""
    Refinement(Tn, CUDAvals, CUDAvecs, Refinement_precision, Multiple_refinement, Max_iter)

Improves the accuracy of the eigenvalues 'CUDAvals' (and vectors 'CUDAvecs') of Tn up to a precision given by Refinement_precision. If Multiple_refinement is true, then will keep refining until it reaches machine precision of the given Refinement_precision, with a maximum of Max_iter iterations.
"""
function Refinement(Tn :: Array, CUDAvals :: Array, CUDAvecs :: Array, Refinement_precision :: Integer, Multiple_refinement :: Bool, Max_iter :: Integer)
    setprecision(BigFloat, Refinement_precision)
    n = length(CUDAvals)
    refinedVals = zeros(BigFloat,n)
    refinedVecs = zeros(BigFloat,n,n)
    Threads.@threads for ii = 1:n
        newVal, newVec = rfeigen(view(Tn, :, :), BigFloat.(view(CUDAvecs, :, ii)), BigFloat(CUDAvals[ii]), BigFloat)
        if Multiple_refinement
            err_norm = abs(dot(view(Tn, 1, :) - newVal, newVec))
            iter = 0
            while err_norm > 100*eps(BigFloat) && iter < Max_iter
                iter += 1
                newVal, newVec = rfeigen(view(Tn, :, :), newVec, newVal, BigFloat)
                err_norm = norm((Tn - newVal*I)*newVec)
            end
        end
        refinedVals[ii] = newVal
        refinedVecs[:,ii] .= newVec
    end
    return refinedVals, refinedVecs
end
