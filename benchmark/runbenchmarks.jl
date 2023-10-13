using PkgBenchmark

config =
    BenchmarkConfig(id=nothing, juliacmd=`julia -O3`, env=Dict("JULIA_NUM_THREADS" => 16))

results = benchmarkpkg("IncrementalInference", config; retune=false)
export_markdown("benchmark/results.md", results)






if false

result = run(SUITE["parametric"])

result = run(SUITE)


foreach(leaves(result)) do bm 
    printstyled("$(bm[1][1]) - $(bm[1][2])\n"; bold=true, reverse=true)
    display(bm[2])
    println()
end


end


if false
## Showing adding variables to fg re-compiles with parametric solve
fg = initfg();
fg.solverParams.graphinit=false;
addVariable!(fg, :x0, Pose2);
addFactor!(fg, [:x0], PriorPose2(MvNormal([0.0,0,0], diagm([0.1,0.1,0.01].^2))));

r = @timed IIF.solveGraphParametric!(fg; init=false, is_sparse=false);

timed = [r];

for i = 1:14
    fr = Symbol("x",i-1)
    to = Symbol("x",i)
    addVariable!(fg, to, Pose2)
    addFactor!(fg, [fr,to], Pose2Pose2(MvNormal([10.0,0,pi/3], diagm([0.5,0.5,0.05].^2))))
    r = @timed IIF.solveGraphParametric!(fg; init=false, is_sparse=false);
    push!(timed, r)
end


addVariable!(fg, :l1, RoME.Point2, tags=[:LANDMARK;]);
addFactor!(fg, [:x0; :l1], Pose2Point2BearingRange(Normal(0.0,0.1), Normal(20.0, 1.0)));
addFactor!(fg, [:x6; :l1], Pose2Point2BearingRange(Normal(0.0,0.1), Normal(20.0, 1.0)));

r = @timed IIF.solveGraphParametric!(fg; init=false, is_sparse=false);
push!(timed, r);

getproperty.(timed, :time)

end




