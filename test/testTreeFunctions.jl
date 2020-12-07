using IncrementalInference
using Test

@testset "Test delete tree clique" begin

fg = generateCanonicalFG_lineStep(3; 
                                   poseEvery=1, 
                                   landmarkEvery=3, 
                                   posePriorsAt=[],
                                   landmarkPriorsAt=[0], 
                                   sightDistance=2,
                                   solverParams=SolverParams(algorithms=[:default, :parametric]))

# getSolverParams(fg).graphinit = false
# getSolverParams(fg).treeinit = true
# getSolverParams(fg).dbg = true
# getSolverParams(fg).useMsgLikelihoods = true
# getSolverParams(fg).dbg = true
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true

smtasks = Task[]
oldtree, smt, hists = solveTree!(fg; smtasks, verbose=true, recordcliqs=ls(fg));

IIF.deleteClique!(oldtree, 1)
# drawTree(oldtree, show=true)

tree, smt, hists = solveTree!(fg, oldtree; smtasks, verbose=true, recordcliqs=ls(fg));

# csmAnimate(tree, hists, frames=1)

end
  