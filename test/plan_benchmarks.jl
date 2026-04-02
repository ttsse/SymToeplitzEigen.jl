using LinearAlgebra
using Printf
using Statistics
using Random
using Dates
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

function _parse_int_list(s::AbstractString)
    vals = Int[]
    for p in split(s, ',')
        t = tryparse(Int, strip(p))
        if !isnothing(t) && t >= 1
            push!(vals, t)
        end
    end
    isempty(vals) && error("No valid integer values were provided.")
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

function _sym_sample_residual(A::AbstractMatrix{Float64}, vals, vecs; sample_count::Integer=3)
    idx = _sample_indices(length(vals), sample_count)
    max_res = 0.0
    for k in idx
        lambda_k = Float64(vals[k])
        xk = Float64.(view(vecs, :, k))
        rk = A * xk .- lambda_k .* xk
        xinf = max(norm(xk, Inf), eps(Float64))
        rel = norm(rk, Inf) / max(abs(lambda_k) * xinf, eps(Float64))
        max_res = max(max_res, rel)
    end
    return max_res
end

function _nonsym_sample_residual(A::AbstractMatrix{Float64}, vals, xr, xl; sample_count::Integer=3)
    idx = _sample_indices(length(vals), sample_count)
    max_res = 0.0
    for k in idx
        lambda_k = vals[k]
        rr = norm((A - lambda_k * I) * xr[:, k], Inf)
        rl = norm((A' - conj(lambda_k) * I) * xl[:, k], Inf)
        max_res = max(max_res, Float64(max(rr, rl)))
    end
    return max_res
end

function _timed_case(runf::Function; repeats::Integer=2, warmups::Integer=1)
    for _ in 1:warmups
        runf()
    end

    times = Float64[]
    bytes = Int[]
    gctimes = Float64[]
    rssdeltas = Int[]
    qualities = Float64[]

    for _ in 1:repeats
        GC.gc()
        rss_before = Sys.maxrss()
        stats = @timed runf()
        rss_after = Sys.maxrss()

        push!(times, stats.time)
        push!(bytes, stats.bytes)
        push!(gctimes, stats.gctime)
        push!(rssdeltas, max(0, Int(rss_after) - Int(rss_before)))

        if stats.value isa Number
            push!(qualities, Float64(stats.value))
        end
    end

    return (
        median_time = median(times),
        min_time = minimum(times),
        median_alloc_bytes = Int(round(median(bytes))),
        median_gc_time = median(gctimes),
        median_peak_rss_delta = Int(round(median(rssdeltas))),
        median_quality = isempty(qualities) ? NaN : median(qualities),
    )
end

function _new_row(; group::String, case::String, threads::Int, n::Int, prec::Int, repeats::Int, stats, speedup::Float64=NaN, efficiency::Float64=NaN, notes::String="")
    return (
        group = group,
        case = case,
        threads = threads,
        n = n,
        prec = prec,
        repeats = repeats,
        median_time = stats.median_time,
        min_time = stats.min_time,
        alloc_mb = stats.median_alloc_bytes / 1024^2,
        peak_rss_mb = stats.median_peak_rss_delta / 1024^2,
        gc_time = stats.median_gc_time,
        quality = stats.median_quality,
        speedup = speedup,
        efficiency = efficiency,
        notes = notes,
    )
end

function _print_rows(rows)
    header = @sprintf(
        "%-24s %-34s %-4s %-6s %-9s %-11s %-11s %-10s %-10s %-12s %-9s %-10s",
        "group", "case", "th", "n", "prec", "med t(s)", "alloc(MB)", "rss(MB)", "gc(s)", "quality", "speedup", "eff"
    )
    println(header)
    println(repeat("-", length(header)))

    for r in rows
        q = isnan(r.quality) ? "-" : @sprintf("%.3e", r.quality)
        sp = isnan(r.speedup) ? "-" : @sprintf("%.2f", r.speedup)
        ef = isnan(r.efficiency) ? "-" : @sprintf("%.2f", r.efficiency)

        println(@sprintf(
            "%-24s %-34s %-4d %-6d %-9d %-11.4f %-11.2f %-10.2f %-10.4f %-12s %-9s %-10s",
            r.group,
            r.case,
            r.threads,
            r.n,
            r.prec,
            r.median_time,
            r.alloc_mb,
            r.peak_rss_mb,
            r.gc_time,
            q,
            sp,
            ef,
        ))
    end
end

function _print_summary(rows)
    function _findrow(case_name::String)
        for r in rows
            if r.case == case_name
                return r
            end
        end
        return nothing
    end

    base = _findrow("baseline_dense_fast")
    opt = _findrow("optimized_auto_fast")

    println()
    println("Summary highlights")
    if !(isnothing(base) || isnothing(opt))
        t_gain = 100 * (base.median_time - opt.median_time) / max(base.median_time, eps(Float64))
        a_gain = 100 * (base.alloc_mb - opt.alloc_mb) / max(base.alloc_mb, eps(Float64))
        println(@sprintf("- Core baseline -> optimized: time %+0.2f%%, alloc %+0.2f%%", t_gain, a_gain))
    end

    println("- Scalability rows report speedup and efficiency against thread-1 optimized run.")
end

function _write_csv(rows, out_path::String)
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, "group,case,threads,n,prec,repeats,median_time_s,min_time_s,alloc_mb,peak_rss_delta_mb,gc_time_s,quality,speedup,efficiency,notes")
        for r in rows
            q = isnan(r.quality) ? "" : @sprintf("%.16e", r.quality)
            sp = isnan(r.speedup) ? "" : @sprintf("%.16e", r.speedup)
            ef = isnan(r.efficiency) ? "" : @sprintf("%.16e", r.efficiency)
            notes = replace(r.notes, "," => ";")
            println(io, @sprintf(
                "%s,%s,%d,%d,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%s,%s,%s,%s",
                r.group,
                r.case,
                r.threads,
                r.n,
                r.prec,
                r.repeats,
                r.median_time,
                r.min_time,
                r.alloc_mb,
                r.peak_rss_mb,
                r.gc_time,
                q,
                sp,
                ef,
                notes,
            ))
        end
    end
end

function _emit_probe_result(stats)
    println(
        "RESULT",
        " time=", @sprintf("%.16e", stats.median_time),
        " min_time=", @sprintf("%.16e", stats.min_time),
        " alloc_bytes=", Int(round(stats.median_alloc_bytes)),
        " peak_rss_delta=", stats.median_peak_rss_delta,
        " gc_time=", @sprintf("%.16e", stats.median_gc_time),
        " quality=", @sprintf("%.16e", stats.median_quality),
    )
end

function _parse_probe_result(out::String)
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
    return (
        median_time = parse(Float64, kv["time"]),
        min_time = parse(Float64, kv["min_time"]),
        median_alloc_bytes = parse(Int, kv["alloc_bytes"]),
        median_peak_rss_delta = parse(Int, kv["peak_rss_delta"]),
        median_gc_time = parse(Float64, kv["gc_time"]),
        median_quality = parse(Float64, kv["quality"]),
    )
end

function _run_scalability_probe(; n::Integer, prec::Integer, max_iter::Integer, repeats::Integer)
    v = [2.0, -1.0]
    A = SymToeplitzEigen.toeplitz(n, v, v)

    runf = () -> begin
        vals, vecs = EigenRef(
            A;
            Refinement_precision = prec,
            Max_iter = max_iter,
            solve_mode = :fast,
            toeplitz_kernel = :auto,
            toeplitz_auto_threshold = 96,
        )
        _sym_sample_residual(A, vals, vecs; sample_count=3)
    end

    return _timed_case(runf; repeats=repeats, warmups=1)
end

function _run_scalability_rows!(rows, threads_list::Vector{Int}; n::Integer, prec::Integer, max_iter::Integer, repeats::Integer)
    tlist = unique(sort(vcat(1, threads_list)))
    probe_rows = NamedTuple[]

    for t in tlist
        cmd = `$(Base.julia_cmd()) --project=$(_PROJECT_ROOT) $(_SCRIPT_PATH) --probe-scalability --n=$n --prec=$prec --max-iter=$max_iter --repeats=$repeats`
        out = read(addenv(cmd, "JULIA_NUM_THREADS" => string(t)), String)
        stats = _parse_probe_result(out)
        push!(probe_rows, (threads=t, stats=stats))
    end

    baseline_t = probe_rows[1].stats.median_time
    for pr in probe_rows
        t = pr.threads
        stats = pr.stats
        speedup = baseline_t / max(stats.median_time, eps(Float64))
        efficiency = speedup / t
        row = _new_row(
            group = "phase2+4+7",
            case = "thread_scaling_auto_fast",
            threads = t,
            n = n,
            prec = prec,
            repeats = repeats,
            stats = stats,
            speedup = speedup,
            efficiency = efficiency,
            notes = "optimized_auto_fast",
        )
        push!(rows, row)
    end
end

function _preset_config(name::AbstractString)
    if lowercase(name) == "quick"
        return (
            core_n = 120,
            policy_n = 40,
            adaptive_n = 40,
            toeplitz_n = 120,
            tao_n = 60,
            nonsym_n = 8,
            scaling_n = 120,
            max_iter = 4,
            repeats = 1,
            cache_calls = 2,
        )
    elseif lowercase(name) == "full"
        return (
            core_n = 5000,
            policy_n = 200,
            adaptive_n = 200,
            toeplitz_n = 5000,
            tao_n = 800,
            nonsym_n = 32,
            scaling_n = 5000,
            max_iter = 40,
            repeats = 2,
            cache_calls = 6,
        )
    end

    throw(ArgumentError("Unsupported preset $(name). Use quick or full."))
end

function run_plan_benchmark_suite(; preset::AbstractString="quick", prec::Integer=256, repeats::Integer=0, max_iter::Integer=0, threads::Vector{Int}=[1, 2, 4], out_path::AbstractString="")
    cfg = _preset_config(preset)
    rep = repeats > 0 ? repeats : cfg.repeats
    iters = max_iter > 0 ? max_iter : cfg.max_iter

    rows = NamedTuple[]

    println("SymToeplitzEigen plan benchmark suite")
    println("preset=", preset, ", prec=", prec, ", repeats=", rep, ", max_iter=", iters)
    println("thread list for scalability=", threads)
    println()

    v = [2.0, -1.0]

    # Group 1: Core evolution across phases 2/3/4.
    A_core = SymToeplitzEigen.toeplitz(cfg.core_n, v, v)
    A_adaptive = SymToeplitzEigen.toeplitz(cfg.adaptive_n, v, v)
    core_cases = [
        ("baseline_dense_fast", A_core, cfg.core_n, (solve_mode=:fast, toeplitz_kernel=:dense, adaptive_precision_escalation=false)),
        ("optimized_auto_fast", A_core, cfg.core_n, (solve_mode=:fast, toeplitz_kernel=:auto, toeplitz_auto_threshold=96, adaptive_precision_escalation=false)),
        ("optimized_adaptive_escalation", A_adaptive, cfg.adaptive_n, (solve_mode=:adaptive, toeplitz_kernel=:auto, toeplitz_auto_threshold=96, adaptive_precision_escalation=true, escalation_precision=2 * prec, escalation_extra_iter=8)),
    ]

    for (name, A_case, n_case, kwargs) in core_cases
        runf = () -> begin
            vals, vecs = EigenRef(A_case; Refinement_precision=prec, Max_iter=iters, kwargs...)
            _sym_sample_residual(A_case, vals, vecs; sample_count=3)
        end
        stats = _timed_case(runf; repeats=rep, warmups=1)
        push!(rows, _new_row(group="phase2+3+4", case=name, threads=Threads.nthreads(), n=n_case, prec=prec, repeats=rep, stats=stats))
    end

    # Group 2: Solve policy comparison (phase 3).
    A_policy = SymToeplitzEigen.toeplitz(cfg.policy_n, v, v)
    policy_cases = [
        ("policy_fast_dense", (solve_mode=:fast, toeplitz_kernel=:dense)),
        ("policy_robust_dense", (solve_mode=:robust, toeplitz_kernel=:dense)),
        ("policy_adaptive_dense", (solve_mode=:adaptive, toeplitz_kernel=:dense, stall_ratio=0.95, stall_iters=2)),
    ]

    for (name, kwargs) in policy_cases
        runf = () -> begin
            vals, vecs = EigenRef(A_policy; Refinement_precision=prec, Max_iter=iters, kwargs...)
            _sym_sample_residual(A_policy, vals, vecs; sample_count=3)
        end
        stats = _timed_case(runf; repeats=rep, warmups=1)
        push!(rows, _new_row(group="phase3", case=name, threads=Threads.nthreads(), n=cfg.policy_n, prec=prec, repeats=rep, stats=stats))
    end

    # Group 3: Toeplitz kernels + cache behavior (phase 4).
    A_toep = SymToeplitzEigen.toeplitz(cfg.toeplitz_n, v, v)
    kernel_cases = [
        ("kernel_dense", (toeplitz_kernel=:dense,)),
        ("kernel_structured", (toeplitz_kernel=:structured,)),
        ("kernel_auto", (toeplitz_kernel=:auto, toeplitz_auto_threshold=96)),
    ]

    for (name, kwargs) in kernel_cases
        runf = () -> begin
            vals, vecs = EigenRef(A_toep; Refinement_precision=prec, Max_iter=iters, solve_mode=:fast, kwargs...)
            _sym_sample_residual(A_toep, vals, vecs; sample_count=3)
        end
        stats = _timed_case(runf; repeats=rep, warmups=1)
        push!(rows, _new_row(group="phase4", case=name, threads=Threads.nthreads(), n=cfg.toeplitz_n, prec=prec, repeats=rep, stats=stats))
    end

    cache_cases = [
        ("cache_off_repeated", false),
        ("cache_on_repeated", true),
    ]

    for (name, use_cache) in cache_cases
        runf = () -> begin
            q = 0.0
            for _ in 1:cfg.cache_calls
                vals, vecs = EigenRef(
                    cfg.toeplitz_n,
                    v;
                    Refinement_precision = prec,
                    Max_iter = max(6, iters ÷ 2),
                    solve_mode = :fast,
                    toeplitz_kernel = :auto,
                    toeplitz_auto_threshold = 96,
                    reuse_toeplitz_cache = use_cache,
                )
                q = _sym_sample_residual(A_toep, vals, vecs; sample_count=2)
            end
            q
        end
        stats = _timed_case(runf; repeats=rep, warmups=1)
        push!(rows, _new_row(group="phase4", case=name, threads=Threads.nthreads(), n=cfg.toeplitz_n, prec=prec, repeats=rep, stats=stats, notes="calls=$(cfg.cache_calls)"))
    end

    # Group 4: Tao optional track overhead (phase 6).
    A_tao = SymToeplitzEigen.toeplitz(cfg.tao_n, v, v)
    tao_cases = [
        ("tao_off", NamedTuple()),
        ("tao_scale_init", (tao_scale_init=true,)),
        ("tao_check", (return_status=true, tao_check=true, tao_max_pairs=10)),
        ("tao_both", (return_status=true, tao_check=true, tao_scale_init=true, tao_max_pairs=10)),
    ]

    for (name, kwargs) in tao_cases
        runf = () -> begin
            if haskey(kwargs, :return_status)
                vals, vecs, _ = EigenRef(A_tao; Refinement_precision=prec, Max_iter=max(8, iters ÷ 2), solve_mode=:fast, toeplitz_kernel=:auto, kwargs...)
                return _sym_sample_residual(A_tao, vals, vecs; sample_count=3)
            else
                vals, vecs = EigenRef(A_tao; Refinement_precision=prec, Max_iter=max(8, iters ÷ 2), solve_mode=:fast, toeplitz_kernel=:auto, kwargs...)
                return _sym_sample_residual(A_tao, vals, vecs; sample_count=3)
            end
        end
        stats = _timed_case(runf; repeats=rep, warmups=1)
        push!(rows, _new_row(group="phase6", case=name, threads=Threads.nthreads(), n=cfg.tao_n, prec=prec, repeats=rep, stats=stats))
    end

    # Group 5: Nonsymmetric path and status overhead (phase 5).
    Random.seed!(12345)
    A_nonsym = randn(cfg.nonsym_n, cfg.nonsym_n)
    nonsym_cases = [
        ("nonsym_no_status", false),
        ("nonsym_with_status", true),
    ]

    for (name, with_status) in nonsym_cases
        runf = () -> begin
            if with_status
                vals, xr, xl, _ = EigenRefNonSym(A_nonsym; Refinement_precision=prec, Max_iter=max(10, iters), return_status=true)
                return _nonsym_sample_residual(A_nonsym, vals, xr, xl; sample_count=3)
            else
                vals, xr, xl = EigenRefNonSym(A_nonsym; Refinement_precision=prec, Max_iter=max(10, iters), return_status=false)
                return _nonsym_sample_residual(A_nonsym, vals, xr, xl; sample_count=3)
            end
        end
        stats = _timed_case(runf; repeats=rep, warmups=1)
        push!(rows, _new_row(group="phase5", case=name, threads=Threads.nthreads(), n=cfg.nonsym_n, prec=prec, repeats=rep, stats=stats))
    end

    # Group 6: Thread scalability for optimized symmetric path (phases 2/4/7).
    _run_scalability_rows!(rows, threads; n=cfg.scaling_n, prec=prec, max_iter=iters, repeats=rep)

    _print_rows(rows)
    _print_summary(rows)

    resolved_out = out_path
    if isempty(resolved_out)
        ts = Dates.format(now(), "yyyymmdd_HHMMSS")
        resolved_out = joinpath(_PROJECT_ROOT, "test", "benchmark_results", "plan_benchmark_" * ts * ".csv")
    end
    _write_csv(rows, resolved_out)

    println()
    println("CSV written to ", resolved_out)

    return rows, resolved_out
end

function _main(args::Vector{String})
    if "--probe-scalability" in args
        n = parse(Int, _arg_value(args, "--n", "1000"))
        prec = parse(Int, _arg_value(args, "--prec", "256"))
        max_iter = parse(Int, _arg_value(args, "--max-iter", "20"))
        repeats = parse(Int, _arg_value(args, "--repeats", "1"))
        stats = _run_scalability_probe(; n=n, prec=prec, max_iter=max_iter, repeats=repeats)
        _emit_probe_result(stats)
        return 0
    end

    preset = _arg_value(args, "--preset", "quick")
    prec = parse(Int, _arg_value(args, "--prec", "256"))
    repeats = parse(Int, _arg_value(args, "--repeats", "0"))
    max_iter = parse(Int, _arg_value(args, "--max-iter", "0"))
    threads = _parse_int_list(_arg_value(args, "--threads", "1,2,4"))
    out_path = _arg_value(args, "--out", "")

    run_plan_benchmark_suite(
        preset=preset,
        prec=prec,
        repeats=repeats,
        max_iter=max_iter,
        threads=threads,
        out_path=out_path,
    )

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_main(ARGS))
end
