using SymToeplitzEigen
using Test
using LinearAlgebra
using SparseArrays
using BandedMatrices
using ToeplitzMatrices

# Check whether the non-allocating helper functions
# are performing as expected
@testset "Helper function check" begin

    #Toeplitz test
    v = [2, -1]
    n = 3
    DT = BigFloat
    A = SymToeplitzEigen.toeplitz(n, DT.(v), DT.(v))
    @test A == [2 -1 0; -1 2 -1; 0 -1 2]

    # normalize_Infnorm! test
    x = DT.(rand(n))
    ss = argmax(abs.(x))
    y = normalize(x, Inf)
    s = SymToeplitzEigen.normalize_Infnorm!(x)
    @test s == ss
    @test x ≈ y

    # my_div_vec! test
    y ./= DT(3)
    SymToeplitzEigen.my_div_vec!(n, x, DT(3))
    @test x ≈ y

    # my_neg_mat_vec_mul! test
    r = DT.(zeros(Int64,n))
    SymToeplitzEigen.my_neg_mat_vec_mul!(n, r, A, x)
    @test r ≈ -A*y

    # my_add_diag_elements! test
    B = A + 1I
    C = A + 0I 
    SymToeplitzEigen.my_add_diag_elements!(n, C, DT(1))
    @test C ≈ B

    # my_vec_setvalue_prom! test
    y2 = Float64.(y)
    SymToeplitzEigen.my_vec_setvalue_prom!(n, x, y2)
    @test norm(x - y, Inf) < 10*eps(Float64) 

    # my_vec_setvalue! test
    SymToeplitzEigen.my_vec_setvalue!(n, x, y)
    @test x ≈ y

    # my_add_vec! test
    a = DT.(rand(n))
    b = DT.(rand(n))
    c = DT.(Float64.(a))
    SymToeplitzEigen.my_add_vec!(n, c, b)
    @test c ≈ a .+ b

    # my_add_scaled_vec! test
    α = DT(0.5)
    d = DT.(Float64.(a))
    SymToeplitzEigen.my_add_scaled_vec!(n, d, α, b)
    @test d ≈ a .+ α .* b

    # my_vec_demote! test
    d64 = zeros(Float64, n)
    SymToeplitzEigen.my_vec_demote!(n, d64, d)
    @test d64 ≈ Float64.(d)

    # my_maxabs test
    @test SymToeplitzEigen.my_maxabs(n, d) ≈ norm(d, Inf)

    # my_mat_setvalue! test
    M1 = DT.(rand(n, n))
    M2 = DT.(zeros(Int64, n, n))
    SymToeplitzEigen.my_mat_setvalue!(n, M2, M1)
    @test M2 ≈ M1

    # my_neg_sym_toeplitz_mat_vec! test
    r2 = DT.(zeros(Int64, n))
    ctoe = DT.([2, -1])
    SymToeplitzEigen.my_neg_sym_toeplitz_mat_vec!(n, r2, ctoe, x)
    @test r2 ≈ -A*x

    # my_fill_sym_toeplitz! test (Float64 and BigFloat)
    T64 = zeros(Float64, n, n)
    SymToeplitzEigen.my_fill_sym_toeplitz!(n, T64, Float64.(ctoe))
    @test T64 ≈ Float64.(A)

    Tbf = DT.(zeros(Int64, n, n))
    SymToeplitzEigen.my_fill_sym_toeplitz!(n, Tbf, ctoe)
    @test Tbf ≈ A
end

# Make sure from a single iteration that we dont go "too far away"
# from low accuracy eigenvalues
@testset "Basic Improvement check" begin
    v = [2, -1]
    n = 100
    T = Float64
    A = SymToeplitzEigen.toeplitz(n, T.(v), T.(v))
    vals, vecs = eigen(A)
    nvals, nvecs = EigenRef(A, Max_iter = 1)
    @test norm(vals - nvals, Inf) < 100*eps(Float64)
    @test norm(abs.(normalize(vecs)) - abs.(normalize(nvecs)), Inf) < 100*eps(Float64) 
end

# Check different low precision still approach same solution
@testset "Low precision check" begin
    v = [2, -1]
    n = 100
    T = Float64
    A = SymToeplitzEigen.toeplitz(n, T.(v), T.(v))
    vals32, vecs32 = EigenRef(n, v, Low_pres_type = Float32)
    vals64, vecs64 = EigenRef(n, v, Low_pres_type = Float64)

    @test norm(maximum([norm((A - vals32[kk]I)*vecs32[:,kk], Inf) for kk in 1:n]) - maximum([norm((A - vals64[kk]I)*vecs64[:,kk], Inf) for kk in 1:n]), Inf) < 1e10*eps(BigFloat) 
end

# Check solve policy modes and status reporting
@testset "Solve policy mode check" begin
    v = [2, -1]
    n = 20
    T = Float64
    A = SymToeplitzEigen.toeplitz(n, T.(v), T.(v))

    valsf, vecsf, statusf = EigenRef(A, solve_mode = :fast, return_status = true, Max_iter = 5)
    @test length(statusf.mode_used) == n
    @test length(statusf.iterations_used) == n
    @test length(statusf.converged) == n
    @test length(statusf.promoted_iteration) == n
    @test length(statusf.precision_bits_used) == n
    @test length(statusf.precision_escalated) == n
    @test all(s -> s in (:fast, :fast_maxiter), statusf.mode_used)
    @test statusf.toeplitz_kernel_used in (:dense, :structured)
    @test all(bits -> bits == 256, statusf.precision_bits_used)
    @test all(statusf.precision_escalated .== false)

    valsr, vecsr, statusr = EigenRef(A, solve_mode = :robust, return_status = true, Max_iter = 5)
    @test all(s -> s in (:robust, :robust_maxiter), statusr.mode_used)

    valsa, vecsa, statusa = EigenRef(A, solve_mode = :adaptive, stall_ratio = 1.0, stall_iters = 1, return_status = true, Max_iter = 5)
    @test all(s -> s in (:adaptive_fast, :adaptive_fast_maxiter, :adaptive_promoted, :adaptive_promoted_maxiter, :adaptive_precision_escalated, :adaptive_precision_escalated_maxiter), statusa.mode_used)
    @test all(it -> it >= 0, statusa.promoted_iteration)

    resf = maximum(norm((A - valsf[kk]I) * vecsf[:, kk], Inf) for kk in 1:n)
    resr = maximum(norm((A - valsr[kk]I) * vecsr[:, kk], Inf) for kk in 1:n)
    resa = maximum(norm((A - valsa[kk]I) * vecsa[:, kk], Inf) for kk in 1:n)
    @test resf < 1e12*eps(BigFloat)
    @test resr < 1e12*eps(BigFloat)
    @test resa < 1e12*eps(BigFloat)

    valsae, vecsae, statusae = EigenRef(
        A,
        solve_mode = :adaptive,
        return_status = true,
        Refinement_precision = 96,
        Max_iter = 1,
        adaptive_precision_escalation = true,
        escalation_precision = 192,
        escalation_extra_iter = 3,
    )
    @test any(statusae.precision_escalated)
    @test all(bits -> bits in (96, 192), statusae.precision_bits_used)
    @test all(s -> s in (:adaptive_fast, :adaptive_fast_maxiter, :adaptive_promoted, :adaptive_promoted_maxiter, :adaptive_precision_escalated, :adaptive_precision_escalated_maxiter), statusae.mode_used)
    resa_e = maximum(norm((A - valsae[kk]I) * vecsae[:, kk], Inf) for kk in 1:n)
    @test resa_e < 1e12*eps(BigFloat)

    @test_throws ArgumentError EigenRef(A, solve_mode = :invalid_mode)
    @test_throws ArgumentError EigenRef(A, solve_mode = :adaptive, adaptive_precision_escalation = true, Refinement_precision = 128, escalation_precision = 64)
    @test_throws ArgumentError EigenRef(A, solve_mode = :adaptive, adaptive_precision_escalation = true, escalation_extra_iter = 0)

    vals_prog, vecs_prog = EigenRef(A, Max_iter = 2, show_progress = true)
    @test length(vals_prog) == n
    @test size(vecs_prog) == size(A)
end

@testset "Toeplitz kernel and cache check" begin
    v = [2, -1]
    n = 30
    T = Float64
    A = SymToeplitzEigen.toeplitz(n, T.(v), T.(v))

    vals_dense, vecs_dense = EigenRef(A, toeplitz_kernel = :dense, Max_iter = 5)
    vals_struct, vecs_struct = EigenRef(A, toeplitz_kernel = :structured, Max_iter = 5)
    vals_auto_struct, vecs_auto_struct = EigenRef(A, toeplitz_kernel = :auto, toeplitz_auto_threshold = 1, Max_iter = 5)
    vals_auto_dense, vecs_auto_dense, status_auto_dense = EigenRef(A, toeplitz_kernel = :auto, toeplitz_auto_threshold = 10_000, return_status = true, Max_iter = 5)

    @test norm(vals_dense - vals_struct, Inf) < 1e10*eps(BigFloat)
    @test norm(abs.(normalize(vecs_dense)) - abs.(normalize(vecs_struct)), Inf) < 1e10*eps(BigFloat)
    @test norm(vals_struct - vals_auto_struct, Inf) < 1e10*eps(BigFloat)
    @test norm(abs.(normalize(vecs_struct)) - abs.(normalize(vecs_auto_struct)), Inf) < 1e10*eps(BigFloat)
    @test norm(vals_dense - vals_auto_dense, Inf) < 1e10*eps(BigFloat)
    @test norm(abs.(normalize(vecs_dense)) - abs.(normalize(vecs_auto_dense)), Inf) < 1e10*eps(BigFloat)
    @test status_auto_dense.toeplitz_kernel_used == :dense

    B = Matrix{Float64}(I, n, n)
    B[diagind(B)] .= range(2.0, 3.0, length=n)
    @inbounds begin
        for kk in 1:(n-1)
            B[kk, kk+1] = 0.1
            B[kk+1, kk] = 0.1
        end
    end
    @test issymmetric(B)
    @test_throws ArgumentError EigenRef(B, toeplitz_kernel = :structured, Max_iter = 1)

    vals_non_toep, vecs_non_toep, status_non_toep = EigenRef(B, toeplitz_kernel = :auto, toeplitz_auto_threshold = 1, return_status = true, Max_iter = 1)
    @test status_non_toep.toeplitz_kernel_used == :dense
    @test maximum(norm((B - vals_non_toep[kk]I) * vecs_non_toep[:, kk], Inf) for kk in 1:n) < 1e-20

    vc = Float64.(v)
    T1 = SymToeplitzEigen._cached_toeplitz(n, vc; reuse=true)
    T2 = SymToeplitzEigen._cached_toeplitz(n, vc; reuse=true)
    @test T1 === T2

    T3 = SymToeplitzEigen._cached_toeplitz(n, vc; reuse=false)
    @test !(T1 === T3)

    @test_throws ArgumentError EigenRef(A, toeplitz_kernel = :auto, toeplitz_auto_threshold = 0)
end

@testset "Sparse and extension matrix support check" begin
    n = 80
    A_sparse = spdiagm(-1 => fill(-1.0, n - 1), 0 => fill(2.0, n), 1 => fill(-1.0, n - 1))
    A_dense = Matrix(A_sparse)

    vals_sp, vecs_sp, status_sp = EigenRef(A_sparse, return_status = true, Max_iter = 4, sparse_solver = :sparse)
    @test length(vals_sp) == n
    @test size(vecs_sp) == (n, n)
    @test status_sp.sparse_solver_used == :sparse
    @test maximum(norm((A_dense - vals_sp[kk]I) * vecs_sp[:, kk], Inf) for kk in 1:n) < 1e4*eps(Float64)

    c = zeros(Float64, n)
    c[1] = 2.0
    c[2] = -1.0
    Tm = SymmetricToeplitz(c)
    vals_toep, vecs_toep = EigenRef(Tm, Max_iter = 4)
    vals_dense_ref, vecs_dense_ref = EigenRef(SymToeplitzEigen.toeplitz(n, [2.0, -1.0], [2.0, -1.0]), Max_iter = 4)
    @test norm(vals_toep - vals_dense_ref, Inf) < 1e10*eps(BigFloat)
    @test norm(abs.(normalize(vecs_toep)) - abs.(normalize(vecs_dense_ref)), Inf) < 1e10*eps(BigFloat)

    Bm = BandedMatrix(A_dense, (1, 1))
    vals_band, vecs_band, status_band = EigenRef(Bm, return_status = true, Max_iter = 4, banded_backend = :sparse)
    @test status_band.sparse_solver_used in (:sparse, :dense)
    @test maximum(norm((A_dense - vals_band[kk]I) * vecs_band[:, kk], Inf) for kk in 1:n) < 1e4*eps(Float64)
end

@testset "Nonsymmetric refinement check" begin
    A = [2.0 1.0 0.0;
         0.0 3.0 1.0;
         0.0 0.0 4.0]

    vals, xr, xl, status = EigenRefNonSym(A, Max_iter = 20, return_status = true)

    @test length(vals) == size(A, 1)
    @test size(xr) == size(A)
    @test size(xl) == size(A)
    @test length(status.mode_used) == size(A, 1)
    @test length(status.converged) == size(A, 1)
    @test length(status.right_residual) == size(A, 1)
    @test length(status.left_residual) == size(A, 1)
    @test length(status.backward_error) == size(A, 1)
    @test length(status.biorthogonality_error) == size(A, 1)
    @test length(status.condition_proxy) == size(A, 1)
    @test length(status.warnings) == size(A, 1)
    @test all(s -> s in (:iterative, :iterative_maxiter, :schur_block_fallback, :factorization_failed, :biorthogonality_failed), status.mode_used)

    max_rr = maximum(norm((A - vals[k]I) * xr[:, k], Inf) for k in 1:size(A, 1))
    max_rl = maximum(norm((A' - conj(vals[k])I) * xl[:, k], Inf) for k in 1:size(A, 1))
    @test max_rr < 1e-20
    @test max_rl < 1e-20

    A2 = [1.0 1e-6;
          0.0 1.0]
    vals2, xr2, xl2, status2 = EigenRefNonSym(A2, Max_iter = 10, return_status = true, cluster_gap_factor = 1e12)
    @test any(status2.clustered)
    @test any(status2.mode_used .== :schur_block_fallback)
    fb_idx = findall(status2.mode_used .== :schur_block_fallback)
    @test !isempty(fb_idx)
    @test all(!isempty(status2.warnings[ii]) for ii in fb_idx)
    @test size(xr2) == size(A2)
    @test size(xl2) == size(A2)

    conv_idx = findall(status.converged)
    @test !isempty(conv_idx)
    @test maximum(status.right_residual[conv_idx]) < 1e-20
    @test maximum(status.left_residual[conv_idx]) < 1e-20
    @test maximum(status.biorthogonality_error[conv_idx]) < 1e-20

    vals_prog, xr_prog, xl_prog = EigenRefNonSym(A, Max_iter = 5, show_progress = true)
    @test length(vals_prog) == size(A, 1)
    @test size(xr_prog) == size(A)
    @test size(xl_prog) == size(A)
end

# Check using both low precision types that we approach the correct eigenvalues
# We check using the Laplace matrix which has known eigenvalues and eigenvectors
@testset "Towards exactness check" begin
    v = [2, -1]
    n = 100
    T = Float64
    A = SymToeplitzEigen.toeplitz(n, T.(v), T.(v))

    nvals, nvecs = EigenRef(A)

    tn = LinRange(BigFloat(pi)/(n+1), n*BigFloat(pi)/(n+1), n)
    tvals = 2 .- 2*cos.(tn) #exact eigenvalues
    tvecs = sin.([1:n;]*tn') #exact eigenvectors

    @test norm(tvals - nvals, Inf) < 1e5*eps(BigFloat)
    @test norm(abs.(normalize(tvecs)) - abs.(normalize(nvecs)), Inf) < 1e5*eps(BigFloat)
end

# Check as we increase precision that error does decrease
@testset "Refinement check for prec $prec" for prec in (128, 256, 512)
    v = [2, -1]
    n = 50
    T = Float64
    A = SymToeplitzEigen.toeplitz(n, T.(v), T.(v))

    nvals, nvecs = EigenRef(A, Refinement_precision = prec)

    tn = LinRange(BigFloat(pi)/(n+1), n*BigFloat(pi)/(n+1), n)
    tvals = 2 .- 2*cos.(tn) #exact eigenvalues
    tvecs = sin.([1:n;]*tn') #exact eigenvectors

    @test norm(tvals - nvals, Inf) < 1e5*eps(BigFloat)
    @test norm(maximum([norm((A - tvals[kk]I)*tvecs[:,kk], Inf) for kk in 1:n]) - maximum([norm((A - nvals[kk]I)*nvecs[:,kk], Inf) for kk in 1:n]), Inf) < 1e5*eps(BigFloat)
end