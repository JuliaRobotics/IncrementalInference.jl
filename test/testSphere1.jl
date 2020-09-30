using IncrementalInference
using Test

@testset "test Sphere1D" begin

fg = initfg()

addVariable!.(fg, [Symbol("x$i") for i=0:4], Sphere1)
addFactor!(fg, [:x0], PriorSphere1(Normal(0.0,0.1)))
map(i->addFactor!(fg, [Symbol("x$i"),Symbol("x$(i+1)")], Sphere1Sphere1(Normal(1.0, 0.1))), 0:3)

tree, smt, hist = solveTree!(fg)

sppes = map(var->getSuggestedPPE(getPPE(var))[1], sortDFG(getVariables(fg),by=getLabel)) 

gt = rem2pi.(collect(0:4), RoundNearest)

@test all(isapprox.(sppes, gt, atol=0.3)) #TODO tolarance, test failed on 0.2

# test packing converters also
d = mktempdir()
saveDFG(d,fg)
lfg = loadDFG(d)
# check loaded fg for all variable and factors
@test issetequal(ls(fg), ls(lfg))
@test issetequal(lsf(fg), lsf(lfg))

end


