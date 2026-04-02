using LinearAlgebra
using Printf
using Statistics
using SymToeplitzEigen

"""
    run_case(n; low=Float64, prec=256, max_iter=100, repeats=3, kernel=:dense)

Runs repeated timed trials for one benchmark case and returns summary metrics.
"""
function run_case(n::Integer; low::Type=Float64, prec::Integer=256, max_iter::Integer=100, repeats::Integer=3, kernel::Symbol=:dense)
    v = [2.0, -1.0]
    A = SymToeplitzEigen.toeplitz(n, low.(v), low.(v))

    # Warm-up compilation and first execution outside measured trials.
    EigenRef(A; Low_pres_type=low, Refinement_precision=prec, Max_iter=max_iter, toeplitz_kernel=kernel)

    times = Float64[]
    bytes = Int[]
    gctimes = Float64[]
    maxres = BigFloat[]

    for _ in 1:repeats
        stats = @timed begin
            vals, vecs = EigenRef(A; Low_pres_type=low, Refinement_precision=prec, Max_iter=max_iter, toeplitz_kernel=kernel)
            maximum(norm((A - vals[k] * I) * vecs[:, k], Inf) for k in 1:n)
        end

        push!(times, stats.time)
        push!(bytes, stats.bytes)
        push!(gctimes, stats.gctime)
        push!(maxres, BigFloat(stats.value))
    end

    return (
        n = n,
        low = low,
        kernel = kernel,
        prec = prec,
        repeats = repeats,
        median_time = median(times),
        min_time = minimum(times),
        median_alloc_bytes = Int(round(median(bytes))),
        median_gc_time = median(gctimes),
        median_max_residual = median(maxres),
    )
end

"""
    run_suite(; ns=(1000, 5000, 10000), precs=(128, 256, 512), kernels=(:dense, :structured, :auto), repeats=3)

Runs a baseline suite and prints a compact table for optimization tracking.
"""
function run_suite(; ns=(1000, 5000, 10000), precs=(128, 256, 512), kernels=(:dense, :structured, :auto), repeats::Integer=3)
    println("SymToeplitzEigen baseline")
    println("Julia threads: ", Threads.nthreads())
    println("Cases: n in ", collect(ns), ", precision in ", collect(precs), ", kernels in ", collect(kernels), ", repeats = ", repeats)
    println()

    header = @sprintf("%-8s %-11s %-9s %-10s %-12s %-12s %-12s %-16s", "n", "kernel", "low", "prec(bits)", "med time(s)", "min time(s)", "alloc(MB)", "med max residual")
    println(header)
    println(repeat("-", length(header)))

    for n in ns
        for kernel in kernels
            for prec in precs
                result = run_case(n; low=Float64, prec=prec, repeats=repeats, kernel=kernel)
                alloc_mb = result.median_alloc_bytes / 1024^2
                row = @sprintf(
                    "%-8d %-11s %-9s %-10d %-12.4f %-12.4f %-12.2f %-16.5e",
                    result.n,
                    string(result.kernel),
                    string(result.low),
                    result.prec,
                    result.median_time,
                    result.min_time,
                    alloc_mb,
                    Float64(result.median_max_residual),
                )
                println(row)
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_suite()
end
