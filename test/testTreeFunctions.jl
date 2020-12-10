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
oldtree, smt, hists = solveTree!(fg; smtasks=smtasks, verbose=false, recordcliqs=ls(fg));

@test IIF.isRoot(oldtree, IIF.CliqueId(1))
@test IIF.isRoot(oldtree, IIF.getClique(oldtree,1))
@test !IIF.isRoot(oldtree, IIF.CliqueId(2))
@test !IIF.isRoot(oldtree, IIF.CliqueId(3))

IIF.deleteClique!(oldtree, IIF.CliqueId(1))

@test IIF.isRoot(oldtree, IIF.CliqueId(2))
@test IIF.isRoot(oldtree, IIF.CliqueId(3))

@test IIF.getClique(oldtree, :x0) == IIF.getClique(oldtree, IIF.CliqueId(3)) == IIF.getClique(oldtree, 1)
@test IIF.getClique(oldtree, :x3) == IIF.getClique(oldtree, IIF.CliqueId(2)) == IIF.getClique(oldtree, 2)
 
# drawTree(oldtree, show=true)

tree, smt, hists = solveTree!(fg, oldtree; smtasks=smtasks, verbose=false, recordcliqs=ls(fg));

csmAnimate(tree, hists, frames=1)

end
  