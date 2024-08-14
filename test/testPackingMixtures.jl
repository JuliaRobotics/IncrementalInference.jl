

using IncrementalInference
using DistributedFactorGraphs
using JSON3


##

@testset "Test packing of mixture of distributions, 1498" begin
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

##

pf = packFactor(getFactor(fg, :x0x1f1))

##

pf_ = JSON3.write(pf)


##

saveDFG("/tmp/caesar/test_mixture.tar.gz", fg)
fg_ = loadDFG("/tmp/caesar/test_mixture.tar.gz")

Base.rm("/tmp/caesar/test_mixture.tar.gz")

##
end