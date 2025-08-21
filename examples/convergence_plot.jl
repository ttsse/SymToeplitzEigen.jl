using SymToeplitzEigen
using LinearAlgebra
using Plots

# Build a small Laplace-type Toeplitz matrix
n = 20
v = [2, -1]
Tn = SymToeplitzEigen.toeplitz(n, Float64.(v), Float64.(v))
vals, vecs = eigen(Tn)

# Run refinement with error tracking
vals_ref, vecs_ref, errors = SymToeplitzEigen.Refinement(Tn, vals, vecs, 256, 20, 1; track_errors=true)

# Pick one eigenvalue’s error history (say first)
err_hist = errors[1]

# Plot
plot(1:length(err_hist), err_hist, yscale=:log10,
     xlabel="Iteration", ylabel="Error norm",
     title="Convergence of refinement", lw=2, marker=:o)


# Save the plot
savefig("convergence.png")