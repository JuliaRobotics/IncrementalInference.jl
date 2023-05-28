# test for variations in DOF solving,  #227 #316 #430

# using Revise

using IncrementalInference
using Test

##

@testset "Test for variations in N while solving" begin

##

fg = generateGraph_CaesarRing1D(graphinit=true)

##

try
  pts_ = approxConv(fg, :x0x1f1, :x1, N=101)
  @test length(pts_) == 101
catch
  # allow one retry, vary rarely has consecutive optimization failure
  pts_ = approxConv(fg, :x0x1f1, :x1, N=101)
  @test length(pts_) == 101
end

##

@error "MUST RESTORE SOLVE WITH DIFFERENT SIZE N, see #1722"
if false
# Change to N=150 AFTER constructing the graph, so solver must update the belief sample values during inference
getSolverParams(fg).N = 150
# getSolverParams(fg).multiproc = false
# getSolverParams(fg).async = false
smtasks = Task[]
tree = solveTree!(fg; smtasks, recordcliqs=ls(fg));

##

pts_ = getBelief(fg, :x1) |> getPoints 
@test length(pts_) == 150
@test length(pts_[1]) == 1

##

# test with change for second solve
getSolverParams(fg).N = 200
tree = solveTree!(fg)

pts_ = getBelief(fg, :x1) |> getPoints
@test length(pts_) == 200
@test length(pts_[1]) == 1

println("test making N smaller than current")
getSolverParams(fg).N = 99
tree = solveTree!(fg)

pts_ = getBelief(fg, :x1) |> getPoints
@warn "removing older solve N size test, likely to be reviewed and updated to new workflow in the future"
@test length(pts_) == 99
@test length(pts_[1]) == 1
end

##

end

#
