# test for variations in DOF solving,  #227 #316 #430

# using Revise

using IncrementalInference
using Test

##

@testset "Test for variations in N while solving" begin

##

fg = generateCanonicalFG_CaesarRing1D(graphinit=true)

pts_ = approxConv(fg, :x0x1f1, :x1, N=101)
if length(pts_) != 101
  @warn "approxConv not adhering to N=101 != $(length(pts_)), see issue #105"
end

# 150 is non-standard
getSolverParams(fg).N = 150
# getSolverParams(fg).multiproc = false
# getSolverParams(fg).async = false
tree, smt, hist = solveTree!(fg)

pts_ = getBelief(fg, :x1) |> getPoints 
@test length(pts_) == 150
@test length(pts_[1]) == 1

# test with change for second solve
getSolverParams(fg).N = 200
tree, smt, hist = solveTree!(fg)

pts_ = getBelief(fg, :x1) |> getPoints
@test length(pts_) == 200
@test length(pts_[1]) == 1


getSolverParams(fg).N = 99
tree, smt, hist = solveTree!(fg)

pts_ = getBelief(fg, :x1) |> getPoints
@warn "removing older solve N size test, likely to be reviewed and updated to new workflow in the future"
@test_broken length(pts_) == 99
@test length(pts_[1]) == 1

##

end

#
