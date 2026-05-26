# test basic forward convolve, see IIF issue #477

# using Revise

using Test
using IncrementalInference
using Statistics
using TensorCast

##

@testset "Test basic convolution result (#477)..." begin
## Start

fg = initfg()

# first numerical values -- samples from the marginal of X0
addVariable!(fg, :x0, ContinuousScalar)
z1 = Normal(0,0.1)
X0 = [rand(z1, 1) for _ in 1:100]
initVariable!(fg, :x0, X0)  # NOTE, manual init without adding a prior to fg

## predict -- project / conv
addVariable!(fg, :x1, ContinuousScalar)

# 0 -> 1 seconds
# make approx function
z2 = Normal(11,1.0) # odo
statemodel = LinearRelative( z2 )
f = addFactor!(fg, [:x0;:x1], statemodel)

##

X1_ = approxConv(fg, getLabel(f), :x1)

@test 10 < Statistics.mean(getindex.(X1_, 1)) < 12
@test 0.5 < Statistics.std(getindex.(X1_, 1)) < 1.5

@error "Add a test to ensure that approxConv does NOT change target variable state.belief values!"

## measure -- product of beliefs, using `ApproxManifoldProducts.jl`

predX1 = HomotopyDensity_legacy(ContinuousScalar(), X1_)
z3 = Normal(9.5,0.75)
measX1 = HomotopyDensity_legacy(ContinuousScalar(), [rand(z3,1) for _ in 1:100])

# do actual product
posterioriX1 = predX1 * measX1
X1 = getPoints(posterioriX1)

## predict, 1->2 seconds

initVariable!(fg, :x2, X1)  # NOTE, manual init without adding a prior to fg
addVariable!(fg, :x2, ContinuousScalar)


z4 = Normal(8,2.0) # odo
statemodel = LinearRelative( z4 )
f = addFactor!(fg, [:x1;:x2], statemodel)
X2_ = approxConv(fg, getLabel(f), :x2)

@cast X2_[i,j] := X2__[j][i]

@test size(X2_) == (1,100)
@test 15 < Statistics.mean(X2_) < 25

##

end


