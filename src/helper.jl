using LinearAlgebra

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



# function toeplitzBanded(n :: Integer, vc :: Array{T, 1}, vr :: Array{T, 1}) where T <: Number
#     Tn = BandedMatrix(Zeros(T,n,n), (length(vc)-1, length(vr)-1))
#     fill!.(view.((Tn,), [diagind(Tn, 1-ii) for ii in eachindex(vc)]), view(vc, :))
#     fill!.(view.((Tn,), [diagind(Tn, jj-1) for jj in Iterators.drop(eachindex(vr), 1)]), @view vr[2:end])
#     return Tn
# end
