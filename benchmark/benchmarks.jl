
using BenchmarkTools

using IncrementalInference
using RoME

##

# Define a parent BenchmarkGroup to contain our SUITE
const SUITE = BenchmarkGroup()

# Add some child groups to our benchmark SUITE.
SUITE["parametric"] = BenchmarkGroup(
    "1-init" => BenchmarkGroup([1, "par"]),
    "2-solve" => BenchmarkGroup([2, "par"]),
    "3-grow" => BenchmarkGroup([3, "par"]),
)
# SUITE["non-parametric"] = BenchmarkGroup(["2-solve"])
# SUITE["construct"] = BenchmarkGroup()

SUITE["parametric"]["1-init"]["hex"] = @benchmarkable(
    IIF.autoinitParametric!(fg);
    samples = 2,
    seconds = 90,
    setup=(println("1-init fg"); fg=generateGraph_Hexagonal(;graphinit=false, landmark=false))
)

SUITE["parametric"]["2-solve"]["hex"] = @benchmarkable(
    IIF.solveGraphParametric!(fg; init=false);
    samples = 2,
    seconds = 90,
    setup=(println("2-fg-1 solve"); fg=generateGraph_Hexagonal(;graphinit=false, landmark=false))
)
    
SUITE["parametric"]["3-grow"]["hex"] = @benchmarkable(
    IIF.solveGraphParametric!(fg; init=false);
    samples = 2,
    seconds = 90,
    setup=(println("3-fg-2 solve"); fg=generateGraph_Hexagonal(;graphinit=false, landmark=true))
)

SUITE["mmisam"] = BenchmarkGroup(
    "1-init" => BenchmarkGroup([1, "non-par"]),
    "2-solve" => BenchmarkGroup([2, "non-par"]),
    "3-grow" => BenchmarkGroup([3, "non-par"]),
)

SUITE["mmisam"]["2-solve"]["hex"] = @benchmarkable(
    solveGraph!(fg);
    samples = 2,
    seconds = 90,
    setup=(println("fg-1 solve"); fg=generateGraph_Hexagonal(;graphinit=true, landmark=false))
)
    
SUITE["mmisam"]["3-grow"]["hex"] = @benchmarkable(
    solveGraph!(fg);
    samples = 2,
    seconds = 90,
    setup=(println("fg-2 solve"); fg=generateGraph_Hexagonal(;graphinit=true, landmark=true))
)

# TODO maintain order (numbered for now), it's a Dict so not guarantied
leaves(SUITE)

# # If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# # the benchmark parameters. Reusing cached parameters is faster and more reliable
# # than re-tuning `SUITE` every time the file is included.
# paramspath = joinpath(dirname(@__FILE__), "params.json")

# if isfile(paramspath)
#     loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals);
# else
#     tune!(SUITE)
#     BenchmarkTools.save(paramspath, BenchmarkTools.params(SUITE));
# end