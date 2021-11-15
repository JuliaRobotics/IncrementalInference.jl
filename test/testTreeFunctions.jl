using IncrementalInference
using Test

@testset "Test delete tree clique" begin

fg = generateGraph_LineStep(3; 
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
oldtree = solveTree!(fg; smtasks=smtasks, verbose=false, recordcliqs=ls(fg));

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

tree = solveTree!(fg, oldtree; smtasks=smtasks, verbose=false, recordcliqs=ls(fg));

# csmAnimate(tree, hists, frames=1)

end


@testset "test tree show and list functions" begin
  
##

fg = generateGraph_Kaess(graphinit=false)

eo = [:l2; :l1; :x1; :x2; :x3]

tree = buildTreeReset!(fg, eo)

show(tree)

show(tree[:x3])
show(tree[:x1])
show(tree[:l2])

rtl = ls(tree)

# test the return is same structure and relatively sorted
@test 3 == length(rtl)
@test rtl[1][1] == 1
@test intersect(rtl[1][2], [:x3; :x2]) |> length == 2
@test rtl[2][1] == 2
@test intersect(rtl[2][2], [:x1; :l1]) |> length == 2
@test rtl[3][1] == 3
@test intersect(rtl[3][2], [:l2;]) |> length == 1

# test individual tree list

# the root clique with two children
rtl1 = ls(tree, 1)
rtl1 = ls(tree, :x3)
@test rtl1.parent |> length == 0
@test rtl1.children |> length == 2
@test intersect( (x->x[1]).(rtl1.children), [2,3]) |> length == 2
@test intersect(rtl1.children[1][2], [:x1; :l1]) |> length == 2
@test intersect(rtl1.children[2][2], [:l2]) |> length == 1


# first leaf clique with a parent
rtl2 = ls(tree, 2)
rtl2 = ls(tree, :x1)
rtl2 = ls(tree, :l1)

@test rtl2.parent |> length == 1
@test rtl2.children |> length == 0
@test intersect( (x->x[1]).(rtl2.parent), [1]) |> length == 1
@test length(rtl2.children) == 0

# second leaf clique with a parent
rtl3 = ls(tree, 3)
rtl3 = ls(tree, :l2)

@test rtl3.parent |> length == 1
@test rtl3.children |> length == 0
@test intersect( (x->x[1]).(rtl3.parent), [1]) |> length == 1
@test length(rtl3.children) == 0


##

end




##