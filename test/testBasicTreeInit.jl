
using Test
using IncrementalInference

@testset "basic per clique stopping criteria" begin

fg = generateCanonicalFG_lineStep(1)
tree, smt, hist = solveTree!(fg, recordcliqs=[:x0;], limititercliqs=[(:x0=>2);])

@test haskey(hist, 1)

@test hist[1] |> length == 7

end



@testset "test endless cycle case, issue #754" begin


fg = generateCanonicalFG_lineStep(5; 
                                  poseEvery=1, 
                                  landmarkEvery=5, 
                                  posePriorsAt=[0,2], 
                                  sightDistance=4,
                                  solverParams=SolverParams(algorithms=[:default, :parametric]))
                                  
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).limititers = 50
smtasks = Task[]
tree, smt, hist = solveTree!(fg; smtasks=smtasks, verbose=true, timeout=50, recordcliqs=ls(fg));

end



@testset "basic test for tree initialization functionality" begin

# small canonical factor graph, without graphinit
fg = generateCanonicalFG_CaesarRing1D(graphinit=false)
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
@show getLogPath(fg)


# # debug only
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true
# getSolverParams(fg).dbg = true

# mkpath(getLogPath(fg))
# verbosefid = open(joinLogPath(fg, "csmVerbose.log"),"w")

tree, smt, hist = solveTree!(fg, timeout=70) # , verbose=true, verbosefid=verbosefid)

# flush(verbosefid)
# close(verbosefid)
# open(joinLogPath(fg, "csmLogicalReconstructMax.log"),"w") do io
#   IIF.reconstructCSMHistoryLogical(getLogPath(fg), fid=io)
# end


end



#
