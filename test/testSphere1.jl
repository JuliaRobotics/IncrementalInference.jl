using IncrementalInference
using Test

##

@testset "test Sphere1D" begin

##

fg = initfg()
getSolverParams(fg).useMsgLikelihoods = true

addVariable!.(fg, [Symbol("x$i") for i=0:4], Sphere1)
addFactor!(fg, [:x0], PriorSphere1(Normal(0.0,0.1)))
map(i->addFactor!(fg, [Symbol("x$i"),Symbol("x$(i+1)")], Sphere1Sphere1(Normal(1.0, 0.1))), 0:3)


solveTree!(fg);

##

sppes = map(var->getPPE(var).suggested[1], sortDFG(getVariables(fg),by=getLabel)) 

gt = rem2pi.(collect(0:4), RoundNearest)

@show sppes
@test all(isapprox.(sppes, gt, atol=0.35))

# test packing converters also
d = "/tmp/caesar/random/testfg"
saveDFG(d,fg)
lfg = loadDFG(d)
Base.rm(d*".tar.gz")
# check loaded fg for all variable and factors
@test issetequal(ls(fg), ls(lfg))
@test issetequal(lsf(fg), lsf(lfg))

##

end


