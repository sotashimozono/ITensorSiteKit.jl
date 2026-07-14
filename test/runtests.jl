ENV["GKSwstype"] = "100"

using ITensorSiteKit, Test

# Canonical universe + completeness guard — the single source of truth, shared
# VERBATIM with the shard planner (test/ci/plan_shards.jl).
include(joinpath(@__DIR__, "ci", "universe.jl"))

# ── Test selection: FILES > SHARD > ALL ──────────────────────────────
#
#   ITENSORSITEKIT_TEST_FILES="base/test_a.jl,base/test_b.jl"
#       explicit list emitted by the LPT shard planner. MUST be a subset of the
#       canonical universe (the planner cannot smuggle in an unglobbed file).
#   ITENSORSITEKIT_TEST_SHARD="k/N"
#       round-robin shard — timing-agnostic; the planner's fallback and a manual knob.
#   neither
#       run everything.
#
# NOTHING is set by default, so a bare `Pkg.test()` — locally, or from any tool that
# has not adopted the sharded workflow — behaves EXACTLY as before: all files. The
# sharding is additive; it never changes the default meaning of the suite.
#
# ITensorSiteKit has no Aqua suite (Aqua is not a test dep), so the planner's `aqua`
# flag has no consumer here — it survives in the emitted plan only because the planner
# is shared VERBATIM across the fleet.
const _files_spec = get(ENV, "ITENSORSITEKIT_TEST_FILES", "")
const _shard_spec = get(ENV, "ITENSORSITEKIT_TEST_SHARD", "")

const _selected, _mode_desc = if !isempty(_files_spec)
    want = [strip(x) for x in split(_files_spec, ",") if !isempty(strip(x))]
    idx = Dict(test_file_key(d, f) => (d, f) for (d, f) in ALL_TEST_FILES)
    sel = Tuple{String,String}[]
    unknown = String[]
    for w in want
        haskey(idx, w) ? push!(sel, idx[w]) : push!(unknown, String(w))
    end
    isempty(unknown) || error(
        "ITENSORSITEKIT_TEST_FILES lists files outside the canonical universe " *
        "(the planner must only emit globbed files): $(unknown)",
    )
    (sel, "FILES (n=$(length(sel)))")
elseif !isempty(_shard_spec)
    parts = split(_shard_spec, "/")
    length(parts) == 2 ||
        error("ITENSORSITEKIT_TEST_SHARD must be \"k/N\"; got $(repr(_shard_spec))")
    k = tryparse(Int, strip(parts[1]))
    n = tryparse(Int, strip(parts[2]))
    (k !== nothing && n !== nothing) || error(
        "ITENSORSITEKIT_TEST_SHARD must be integer \"k/N\"; got $(repr(_shard_spec))",
    )
    (1 <= k <= n) || error("ITENSORSITEKIT_TEST_SHARD k/N needs 1 ≤ k ≤ N; got $k/$n")
    n <= length(ALL_TEST_FILES) || error(
        "ITENSORSITEKIT_TEST_SHARD N=$n exceeds the $(length(ALL_TEST_FILES))-file " *
        "suite; shards $(length(ALL_TEST_FILES) + 1)..$n would run zero tests — lower N.",
    )
    sel = [tf for (i, tf) in enumerate(ALL_TEST_FILES) if ((i - 1) % n) + 1 == k]
    (sel, "SHARD $k/$n")
else
    (ALL_TEST_FILES, "ALL")
end

println(
    "Test selection: $(_mode_desc) → $(length(_selected))/$(length(ALL_TEST_FILES)) files"
)

# Per-file wall-clock, emitted so the NEXT run's planner can LPT bin-pack instead of
# falling back to round-robin. `julia-actions/julia-runtest` SANDBOXES Pkg.test into a
# temp copy, so `@__DIR__` points somewhere `upload-artifact` will never look — CI pins
# the output dir into the workspace via ITENSORSITEKIT_CIOUT_DIR.
const _emit = get(ENV, "ITENSORSITEKIT_EMIT", "0") == "1"
const _ciout = get(ENV, "ITENSORSITEKIT_CIOUT_DIR", joinpath(@__DIR__, ".ci-out"))
const _sid = get(ENV, "ITENSORSITEKIT_SHARD_ID", "local")
_emit && mkpath(_ciout)

const FIG_BASE = joinpath(pkgdir(ITensorSiteKit), "docs", "src", "assets")
const PATHS = Dict()
mkpath.(values(PATHS))

@testset "tests" begin
    test_args = copy(ARGS)
    println("Passed arguments ARGS = $(test_args) to tests.")

    timings = Tuple{String,Float64}[]
    @time for (d, f) in _selected
        filepath = joinpath(@__DIR__, d, f)
        key = test_file_key(d, f)
        @testset "$f" begin
            t = @elapsed begin
                println("  Including $(filepath)")
                include(filepath)
            end
            println("  [time] $(key): $(round(t; digits=2))s")
            push!(timings, (key, t))
        end
    end

    if _emit
        open(joinpath(_ciout, "timings-$(_sid).tsv"), "w") do io
            for (k, t) in timings
                println(io, k, '\t', round(t; digits=3))
            end
        end
        println("Emitted $(length(timings)) timing rows -> $(_ciout)/timings-$(_sid).tsv")
    end
end
