
using Test
using IncrementalInference
# using Gadfly

##
@testset "basic per clique stopping criteria" begin
##

fg = generateGraph_LineStep(1)
smtasks = Task[]
tree = solveTree!(fg, smtasks=smtasks, recordcliqs=[:x0;], limititercliqs=[(:x0=>2);])
hist = fetchCliqHistoryAll!(smtasks)

@test haskey(hist, 1)

@test hist[1] |> length == 2

#normal solve should have 11 states, update when more are added.
fg = generateGraph_LineStep(1)
smtasks = Task[]
tree = solveTree!(fg, smtasks=smtasks, recordcliqs=[:x0;]);
hist = fetchCliqHistoryAll!(smtasks)

@test haskey(hist, 1)

@test hist[1] |> length == 12

end



@testset "test endless cycle case, issue #754" begin


fg = generateGraph_LineStep(5; 
                            poseEvery=1, 
                            landmarkEvery=5, 
                            posePriorsAt=[0,2], 
                            sightDistance=4,
                            solverParams=SolverParams(algorithms=[:default, :parametric]))
                            
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).limititers = 50
smtasks = Task[]
tree = solveTree!(fg; smtasks=smtasks, verbose=true, timeout=50, recordcliqs=ls(fg));

end



@testset "basic test for tree initialization functionality" begin

# small canonical factor graph, without graphinit
fg = generateGraph_CaesarRing1D(graphinit=false)
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
@show getLogPath(fg)


# # debug only
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true
# getSolverParams(fg).dbg = true

# mkpath(getLogPath(fg))
# verbosefid = open(joinLogPath(fg, "csmVerbose.log"),"w")

tree = solveTree!(fg, timeout=70) # , verbose=true, verbosefid=verbosefid)

# flush(verbosefid)
# close(verbosefid)
# open(joinLogPath(fg, "csmLogicalReconstructMax.log"),"w") do io
#   IIF.reconstructCSMHistoryLogical(getLogPath(fg), fid=io)
# end


end

@testset "basic tree initialization limittreeinit_iters" begin

# part of fg that can init
fg = generateGraph_LineStep(3; poseEvery=1)

good_vars = sortDFG(ls(fg))
# fg.solverParams.showtree = true
# fg.solverParams.drawtree = true

# part of fg that cannot init
addVariable!(fg, :s0, ContinuousScalar)
addVariable!(fg, :s1, ContinuousScalar)
addVariable!(fg, :s2, ContinuousScalar)

addFactor!(fg, [:s0;:s1], LinearRelative(Normal()))
addFactor!(fg, [:s1;:s2], LinearRelative(Normal()))

smtasks = Task[]

tree = solveTree!(fg; smtasks=smtasks, verbose=true)


for var in good_vars
    sppe = calcMeanMaxSuggested(fg,var).suggested
    println("Testing ", var,": ", sppe)
    @test isapprox(sppe[1], parse(Int,string(var)[end]), atol=0.15)
end


@error "Restore test on GadflyExt.spyCliqMat"
# pl = spyCliqMat(getClique(tree,1));

##
end

#
