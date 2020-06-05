# test for variations in DOF solving,  #227 #316 #430

# using Revise

using IncrementalInference
using Test


@testset "Test for variations in N while solving" begin

fg = generateCanonicalFG_CaesarRing1D(graphinit=true)


pts = approxConv(fg, :x0x1f1, :x1, N=101)
if size(pts,2) != 101
  @warn "approxConv not adhering to N=101 != $(size(pts,2)), see issue #105"
end

# 150 is non-standard
getSolverParams(fg).N = 150
# getSolverParams(fg).multiproc = false
# getSolverParams(fg).async = false
tree, smt, hist = solveTree!(fg)

@test getKDE(fg, :x1) |> getPoints |> size == (1,150)

# test with change for second solve
getSolverParams(fg).N = 200
tree, smt, hist = solveTree!(fg)

@test getKDE(fg, :x1) |> getPoints |> size == (1,200)


getSolverParams(fg).N = 99
tree, smt, hist = solveTree!(fg)

@test getKDE(fg, :x1) |> getPoints |> size == (1,99)

end

#
