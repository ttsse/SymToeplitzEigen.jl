"""
    _progress_init(total; show=false, label="Progress")

Internal progress helper.
If `ProgressMeter` is loaded in `Main`, it is used; otherwise a lightweight
text progress bar is printed to `stderr`.
"""
function _progress_init(total::Integer; show::Bool=false, label::AbstractString="Progress")
    if !show || total <= 0
        return nothing
    end

    total_i = Int(total)
    lock = ReentrantLock()
    counter = Threads.Atomic{Int}(0)
    step = max(1, cld(total_i, 50))
    start_time = time()

    pm_mod = nothing
    pm_obj = nothing

    if isdefined(Main, :ProgressMeter)
        m = getfield(Main, :ProgressMeter)
        if m isa Module && isdefined(m, :Progress)
            try
                pm_obj = getfield(m, :Progress)(total_i; desc=string(label, ": "), dt=0.2)
                pm_mod = m
            catch
                pm_mod = nothing
                pm_obj = nothing
            end
        end
    end

    return (
        total = total_i,
        label = String(label),
        lock = lock,
        counter = counter,
        step = step,
        start_time = start_time,
        pm_mod = pm_mod,
        pm_obj = pm_obj,
    )
end

"""
    _progress_update!(state)

Internal progress update called once per completed unit of work.
"""
function _progress_update!(state)
    if isnothing(state)
        return
    end

    done = Threads.atomic_add!(state.counter, 1) + 1

    if !isnothing(state.pm_obj)
        lock(state.lock) do
            getfield(state.pm_mod, :next!)(state.pm_obj)
        end
        return
    end

    if done == state.total || done % state.step == 0
        frac = done / state.total
        filled = clamp(round(Int, 24 * frac), 0, 24)
        bar = "[" * repeat("=", filled) * repeat(" ", 24 - filled) * "]"
        pct = round(100 * frac; digits=1)
        elapsed = round(time() - state.start_time; digits=1)

        lock(state.lock) do
            print(stderr, "\r", state.label, " ", bar, " ", done, "/", state.total, " ", pct, "% ", elapsed, "s")
            if done == state.total
                print(stderr, "\n")
            end
            flush(stderr)
        end
    end
end

"""
    _progress_finish!(state)

Internal progress finalization.
"""
function _progress_finish!(state)
    if isnothing(state)
        return
    end

    if !isnothing(state.pm_obj)
        lock(state.lock) do
            if isdefined(state.pm_mod, :finish!)
                getfield(state.pm_mod, :finish!)(state.pm_obj)
            else
                print(stderr, "\n")
                flush(stderr)
            end
        end
    end
end
