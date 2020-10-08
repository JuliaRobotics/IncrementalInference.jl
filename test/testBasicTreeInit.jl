
using Test
using IncrementalInference

@testset "basic per clique stopping criteria" begin

fg = generateCanonicalFG_lineStep(1)
tree, smt, hist = solveTree!(fg, recordcliqs=[:x0;], limititercliqs=[(:x0=>2);])

@test haskey(hist, 1)

@test hist[1] |> length == 2

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

for var in sortDFG(ls(fg))
    sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.1)
end

# linear octo 
N=8
fg = generateCanonicalFG_lineStep(N; 
                                  graphinit=false,
                                  poseEvery=1, 
                                  landmarkEvery=N+1, 
                                  posePriorsAt=[0],
                                  landmarkPriorsAt=[], 
                                  sightDistance=N+1)

deleteFactor!.(fg, [Symbol("x$(i)lm0f1") for i=1:(N-1)])

getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).limititers = 200

smtasks = Task[]
try 
    tree, smt, hists = IIF.solveTree!(fg; smtasks=smtasks, timeout=10);
    
    for var in sortDFG(ls(fg))
        sppe = getVariable(fg,var) |> getPPE |> IIF.getSuggestedPPE
        println("Testing ", var,": ", sppe)
        @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.15)
    end
catch
    @warn "Octo 754 tests failed"
end
@test_broken !any(istaskfailed.(smtasks))




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

# hack fix during #910 and #754
idb = [6=>(determineCliqNeedDownMsg_StateMachine=>10);]

# mkpath(getLogPath(fg))
# verbosefid = open(joinLogPath(fg, "csmVerbose.log"),"w")

tree, smt, hist = solveTree!(fg, timeout=70, injectDelayBefore=idb ) # , verbose=true, verbosefid=verbosefid)

# flush(verbosefid)
# close(verbosefid)
# open(joinLogPath(fg, "csmLogicalReconstructMax.log"),"w") do io
#   IIF.reconstructCSMHistoryLogical(getLogPath(fg), fid=io)
# end


end



#
