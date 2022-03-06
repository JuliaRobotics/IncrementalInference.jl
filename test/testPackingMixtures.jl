

using IncrementalInference
using DistributedFactorGraphs
using JSON


##

@testset "Test packing of mixture of distributions" begin
##

# Start with an empty factor graph
fg = initfg()

# add the first node
addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
mmo = Mixture(LinearRelative, 
              (hypo1=Rayleigh(3), hypo2=Uniform(30,55)), 
              [0.4; 0.6])
addFactor!(fg, [:x0, :x1], mmo)
JSON.json(packFactor(fg, getFactor(fg, :x0x1f1)))

##
end