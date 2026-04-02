using LinearAlgebra
using Printf
using Statistics
using SymToeplitzEigen

const _SCRIPT_PATH = @__FILE__
const _PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

function _arg_value(args::Vector{String}, key::String, default::String)
    prefix = key * "="
    for a in args
        if startswith(a, prefix)
            return split(a, "=", limit=2)[2]
        end
    end
    return default
end

function _parse_threads(s::AbstractString)
    parts = split(s, ',')
    vals = Int[]
    for p in parts
        t = tryparse(Int, strip(p))
        if !isnothing(t) && t >= 1
            push!(vals, t)
        end
    end
    isempty(vals) && error("No valid thread counts were provided.")
    return unique(sort(vals))
end

function _sample_indices(n::Integer, count::Integer)
    c = max(1, min(Int(count), Int(n)))
    if c == Int(n)
        return collect(1:Int(n))
    end
    idx = round.(Int, range(1, Int(n), length=c))
    return unique(clamp.(idx, 1, Int(n)))
end

function _sym_kwargs(config::Symbol)
    if config == :baseline
        return (
            solve_mode = :fast,
            toeplitz_kernel = :dense,
            adaptive_precision_escalation = false,
        )
    elseif config == :optimized
        return (
            solve_mode = :fast,
            toeplitz_kernel = :auto,
            toeplitz_auto_threshold = 96,
            adaptive_precision_escalation = false,
        )
    end
    throw(ArgumentError("Unknown symmetric probe config: $(config)"))
end

function _sym_sample_quality(A::AbstractMatrix{Float64}, vals, vecs; sample_count::Integer=3)
    idx = _sample_indices(length(vals), sample_count)
    Ainf = opnorm(A, Inf)

    max_res = 0.0
    max_back = 0.0

    for kk in idx
        λk = Float64(vals[kk])
        xk = Float64.(view(vecs, :, kk))
        rk = A * xk .- λk .* xk
        rinf = norm(rk, Inf)
        xinf = max(norm(xk, Inf), eps(Float64))

        rel_res = rinf / max(abs(λk) * xinf, eps(Float64))
        back_err = rinf / max((Ainf + abs(λk)) * xinf, eps(Float64))

        max_res = max(max_res, rel_res)
        max_back = max(max_back, back_err)
    end

    return max_res, max_back
end

function _sym_probe(; config::Symbol=:baseline, n::Integer=5000, prec::Integer=256, max_iter::Integer=40, sample_count::Integer=3)
    v = [2.0, -1.0]
    A = SymToeplitzEigen.toeplitz(n, v, v)
    kwargs = _sym_kwargs(config)

    # Warm-up this configuration so compile time is not included in gate timing.
    EigenRef(A; Refinement_precision=prec, Max_iter=max_iter, return_status=true, kwargs...)
    GC.gc()

    rss_before = Sys.maxrss()
    stats = @timed EigenRef(A; Refinement_precision=prec, Max_iter=max_iter, return_status=true, kwargs...)
    vals, vecs, status = stats.value
    rss_after = Sys.maxrss()

    max_residual, max_backward = _sym_sample_quality(A, vals, vecs; sample_count=sample_count)

    mid = max(1, length(vals) ÷ 2)
    sample1 = Float64(vals[1])
    sample2 = Float64(vals[mid])
    sample3 = Float64(vals[end])

    return (
        config = config,
        n = Int(n),
        prec = Int(prec),
        max_iter = Int(max_iter),
        time = Float64(stats.time),
        alloc_bytes = Int(stats.bytes),
        gc_time = Float64(stats.gctime),
        peak_rss_delta = max(0, Int(rss_after) - Int(rss_before)),
        max_residual = Float64(max_residual),
        max_backward = Float64(max_backward),
        mean_iter = Float64(mean(status.iterations_used)),
        converged_fraction = Float64(count(status.converged) / length(status.converged)),
        toeplitz_kernel_used = String(status.toeplitz_kernel_used),
        sample1 = sample1,
        sample2 = sample2,
        sample3 = sample3,
    )
end

function _emit_sym_probe_result(r)
    println(
        "RESULT",
        " type=sym",
        " config=", String(r.config),
        " n=", r.n,
        " prec=", r.prec,
        " max_iter=", r.max_iter,
        " time=", @sprintf("%.16e", r.time),
        " alloc_bytes=", r.alloc_bytes,
        " gc_time=", @sprintf("%.16e", r.gc_time),
        " peak_rss_delta=", r.peak_rss_delta,
        " max_residual=", @sprintf("%.16e", r.max_residual),
        " max_backward=", @sprintf("%.16e", r.max_backward),
        " mean_iter=", @sprintf("%.16e", r.mean_iter),
        " converged_fraction=", @sprintf("%.16e", r.converged_fraction),
        " kernel=", r.toeplitz_kernel_used,
        " sample1=", @sprintf("%.16e", r.sample1),
        " sample2=", @sprintf("%.16e", r.sample2),
        " sample3=", @sprintf("%.16e", r.sample3),
    )
end

function _parse_result_line(out::String)
    line = nothing
    for l in split(out, '\n')
        sl = strip(l)
        if startswith(sl, "RESULT ")
            line = sl
        end
    end
    isnothing(line) && error("Could not parse probe output.\n" * out)

    kv = Dict{String, String}()
    for token in split(line)[2:end]
        parts = split(token, "=", limit=2)
        if length(parts) == 2
            kv[parts[1]] = parts[2]
        end
    end
    return kv
end

function _run_sym_probe_subprocess(config::Symbol; n::Integer, prec::Integer, max_iter::Integer, sample_count::Integer, threads::Integer)
    cmd = `$(Base.julia_cmd()) --project=$(_PROJECT_ROOT) $(_SCRIPT_PATH) --probe-sym --config=$(String(config)) --n=$n --prec=$prec --max-iter=$max_iter --samples=$sample_count`
    out = read(addenv(cmd, "JULIA_NUM_THREADS" => string(threads)), String)
    kv = _parse_result_line(out)

    return (
        config = Symbol(kv["config"]),
        n = parse(Int, kv["n"]),
        prec = parse(Int, kv["prec"]),
        max_iter = parse(Int, kv["max_iter"]),
        time = parse(Float64, kv["time"]),
        alloc_bytes = parse(Int, kv["alloc_bytes"]),
        gc_time = parse(Float64, kv["gc_time"]),
        peak_rss_delta = parse(Int, kv["peak_rss_delta"]),
        max_residual = parse(Float64, kv["max_residual"]),
        max_backward = parse(Float64, kv["max_backward"]),
        mean_iter = parse(Float64, kv["mean_iter"]),
        converged_fraction = parse(Float64, kv["converged_fraction"]),
        kernel = kv["kernel"],
        sample1 = parse(Float64, kv["sample1"]),
        sample2 = parse(Float64, kv["sample2"]),
        sample3 = parse(Float64, kv["sample3"]),
        threads = Int(threads),
    )
end

function _run_perf_memory_gate(; n::Integer=5000, prec::Integer=256, max_iter::Integer=40, sample_count::Integer=3, threads::Integer=1)
    baseline = _run_sym_probe_subprocess(:baseline; n=n, prec=prec, max_iter=max_iter, sample_count=sample_count, threads=threads)
    optimized = _run_sym_probe_subprocess(:optimized; n=n, prec=prec, max_iter=max_iter, sample_count=sample_count, threads=threads)

    runtime_ok = optimized.time < baseline.time
    peak_mem_ok = optimized.peak_rss_delta < baseline.peak_rss_delta

    return (
        passed = runtime_ok && peak_mem_ok,
        runtime_ok = runtime_ok,
        peak_mem_ok = peak_mem_ok,
        baseline = baseline,
        optimized = optimized,
    )
end

function _run_sym_quality_gate(; n::Integer=1200, prec::Integer=256, max_iter::Integer=30, sample_count::Integer=3, threads=(1, 2, 4), residual_tol::Float64=1e-10, backward_tol::Float64=1e-10, reproducibility_tol::Float64=1e-11, iter_spread_tol::Float64=1.20)
    probes = [_run_sym_probe_subprocess(:optimized; n=n, prec=prec, max_iter=max_iter, sample_count=sample_count, threads=t) for t in threads]

    residual_ok = all(p.max_residual <= residual_tol for p in probes)
    backward_ok = all(p.max_backward <= backward_tol for p in probes)

    ref = probes[1]
    diffs = Float64[]
    for p in probes
        d = maximum(abs.([
            p.sample1 - ref.sample1,
            p.sample2 - ref.sample2,
            p.sample3 - ref.sample3,
        ]))
        push!(diffs, d)
    end
    reproducibility_ok = all(d <= reproducibility_tol for d in diffs)

    iter_vals = [p.mean_iter for p in probes]
    iter_spread = maximum(iter_vals) / max(minimum(iter_vals), eps(Float64))
    iter_consistent_ok = iter_spread <= iter_spread_tol

    return (
        passed = residual_ok && backward_ok && reproducibility_ok && iter_consistent_ok,
        residual_ok = residual_ok,
        backward_ok = backward_ok,
        reproducibility_ok = reproducibility_ok,
        iter_consistent_ok = iter_consistent_ok,
        iter_spread = iter_spread,
        reproducibility_diffs = diffs,
        probes = probes,
    )
end

function _run_nonsymmetric_gate(; prec::Integer=256, max_iter::Integer=20, residual_tol::Float64=1e-20, biorth_tol::Float64=1e-20)
    A = [2.0 1.0 0.0;
         0.0 3.0 1.0;
         0.0 0.0 4.0]

    vals, xr, xl, status = EigenRefNonSym(A; Refinement_precision=prec, Max_iter=max_iter, return_status=true)

    conv_idx = findall(status.converged)
    conv_nonempty = !isempty(conv_idx)
    right_ok = conv_nonempty && maximum(status.right_residual[conv_idx]) <= residual_tol
    left_ok = conv_nonempty && maximum(status.left_residual[conv_idx]) <= residual_tol
    biorth_ok = conv_nonempty && maximum(status.biorthogonality_error[conv_idx]) <= biorth_tol
    cond_field_ok = length(status.condition_proxy) == size(A, 1)
    warning_field_ok = length(status.warnings) == size(A, 1)

    # Clustered scenario should produce explicit warning metadata.
    A2 = [1.0 1e-6;
          0.0 1.0]
    _, _, _, status2 = EigenRefNonSym(A2; Refinement_precision=prec, Max_iter=max_iter, return_status=true, cluster_gap_factor=1e12)
    fb_idx = findall(status2.mode_used .== :schur_block_fallback)
    warning_reporting_ok = !isempty(fb_idx) && all(!isempty(status2.warnings[ii]) for ii in fb_idx)

    # Secondary external check for right/left residuals on the deterministic case.
    max_rr_external = maximum(norm((A - vals[k] * I) * xr[:, k], Inf) for k in 1:size(A, 1))
    max_rl_external = maximum(norm((A' - conj(vals[k]) * I) * xl[:, k], Inf) for k in 1:size(A, 1))
    external_ok = max_rr_external <= residual_tol && max_rl_external <= residual_tol

    return (
        passed = right_ok && left_ok && biorth_ok && cond_field_ok && warning_field_ok && warning_reporting_ok && external_ok,
        right_ok = right_ok,
        left_ok = left_ok,
        biorth_ok = biorth_ok,
        cond_field_ok = cond_field_ok,
        warning_field_ok = warning_field_ok,
        warning_reporting_ok = warning_reporting_ok,
        external_ok = external_ok,
        max_rr_external = max_rr_external,
        max_rl_external = max_rl_external,
    )
end

function run_phase7_release_gates(; n_perf::Integer=5000, n_quality::Integer=1200, prec::Integer=256, max_iter::Integer=40, sample_count::Integer=3, perf_threads::Integer=1, threads=(1, 2, 4))
    perf_gate = _run_perf_memory_gate(; n=n_perf, prec=prec, max_iter=max_iter, sample_count=sample_count, threads=perf_threads)
    quality_gate = _run_sym_quality_gate(; n=n_quality, prec=prec, max_iter=max(10, max_iter ÷ 2), sample_count=sample_count, threads=threads)
    nonsym_gate = _run_nonsymmetric_gate(; prec=prec, max_iter=max(10, max_iter ÷ 2))

    return (
        passed = perf_gate.passed && quality_gate.passed && nonsym_gate.passed,
        perf_gate = perf_gate,
        quality_gate = quality_gate,
        nonsym_gate = nonsym_gate,
    )
end

function _print_gate_report(report)
    println("Phase 7 Release Gates")
    println("- Overall: ", report.passed ? "PASS" : "FAIL")

    pg = report.perf_gate
    println("- Gate 7.1 runtime+peak-memory: ", pg.passed ? "PASS" : "FAIL")
    println(@sprintf("  baseline  time=%.4fs  peak_rss_delta=%.2fMB  alloc=%.2fMB", pg.baseline.time, pg.baseline.peak_rss_delta / 1024^2, pg.baseline.alloc_bytes / 1024^2))
    println(@sprintf("  optimized time=%.4fs  peak_rss_delta=%.2fMB  alloc=%.2fMB", pg.optimized.time, pg.optimized.peak_rss_delta / 1024^2, pg.optimized.alloc_bytes / 1024^2))

    qg = report.quality_gate
    println("- Gate 7.2 numerical quality+thread reproducibility: ", qg.passed ? "PASS" : "FAIL")
    println("  residual_ok=", qg.residual_ok, ", backward_ok=", qg.backward_ok, ", reproducibility_ok=", qg.reproducibility_ok, ", iter_consistent_ok=", qg.iter_consistent_ok)
    println(@sprintf("  iteration_spread=%.4f", qg.iter_spread))

    ng = report.nonsym_gate
    println("- Gate 7.3 nonsymmetric residual+biorthogonality+warnings: ", ng.passed ? "PASS" : "FAIL")
    println("  right_ok=", ng.right_ok, ", left_ok=", ng.left_ok, ", biorth_ok=", ng.biorth_ok, ", cond_field_ok=", ng.cond_field_ok, ", warning_reporting_ok=", ng.warning_reporting_ok)
    println(@sprintf("  external_max_rr=%.3e  external_max_rl=%.3e", Float64(ng.max_rr_external), Float64(ng.max_rl_external)))
end

function _main(args::Vector{String})
    if "--probe-sym" in args
        config = Symbol(_arg_value(args, "--config", "baseline"))
        n = parse(Int, _arg_value(args, "--n", "5000"))
        prec = parse(Int, _arg_value(args, "--prec", "256"))
        max_iter = parse(Int, _arg_value(args, "--max-iter", "40"))
        sample_count = parse(Int, _arg_value(args, "--samples", "3"))
        result = _sym_probe(; config=config, n=n, prec=prec, max_iter=max_iter, sample_count=sample_count)
        _emit_sym_probe_result(result)
        return 0
    end

    n_perf = parse(Int, _arg_value(args, "--n-perf", "5000"))
    n_quality = parse(Int, _arg_value(args, "--n-quality", "1200"))
    prec = parse(Int, _arg_value(args, "--prec", "256"))
    max_iter = parse(Int, _arg_value(args, "--max-iter", "40"))
    sample_count = parse(Int, _arg_value(args, "--samples", "3"))
    perf_threads = parse(Int, _arg_value(args, "--perf-threads", "1"))
    threads = _parse_threads(_arg_value(args, "--threads", "1,2,4"))

    report = run_phase7_release_gates(
        n_perf=n_perf,
        n_quality=n_quality,
        prec=prec,
        max_iter=max_iter,
        sample_count=sample_count,
        perf_threads=perf_threads,
        threads=Tuple(threads),
    )
    _print_gate_report(report)

    return report.passed ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_main(ARGS))
end
