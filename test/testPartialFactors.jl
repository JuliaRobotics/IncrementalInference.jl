# test aspects specific to partial factors

using IncrementalInference
using Test

@testset "testing basic partial factor functions..." begin

fg = initfg()
addVariable!(fg, :x0, ContinuousScalar)

addFactor!(fg, [:x0], Prior(Normal()) )

fc = getFactor(fg, :x0f1)
@test !isPartial(fc)


fg = initfg()

addVariable!(fg, :x1, ContinuousEuclid{2})

addFactor!(fg, [:x1], PartialPrior(ContinuousEuclid{2},Normal(), (1,)) )

@test isPartial(getFactor(fg, :x1f1))

end
