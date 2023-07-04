using SymToeplitzEigen
using Test
using LinearAlgebra

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
    @test norm(maximum([norm((A - tvals[kk]I)*tvecs[:,kk], Inf) for kk in n]) - maximum([norm((A - nvals[kk]I)*nvecs[:,kk], Inf) for kk in n]), Inf) < 1e5*eps(BigFloat)
end