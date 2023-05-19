# test basic forward convolve, see IIF issue #477

# using Revise

using Test
using IncrementalInference
using Statistics
using TensorCast

##

@testset "Test basic convolution result (#477)..." begin

##

function forwardConvolve(X0::AbstractVector{P}, model) where P

  fg = initfg()

  addVariable!(fg, :x0, ContinuousScalar)
  initVariable!(fg, :x0, X0)
  addVariable!(fg, :x1, ContinuousScalar)

  addFactor!(fg, [:x0;:x1], model)

  ## TODO -- dont use name here, add API to just use z2 here
  return approxConv(fg, :x0x1f1, :x1)
end


## Start

# first numerical values -- samples from the marginal of X0
z1 = Normal(0,0.1)
X0 = [rand(z1, 1) for _ in 1:100]


## predict -- project / conv
# 0 -> 1 seconds
# make approx function
z2 = Normal(11,1.0) # odo
statemodel = LinearRelative( z2 )
X1_ = forwardConvolve(X0, statemodel)


## measure -- product of beliefs, using `ApproxManifoldProducts.jl`

predX1 = manikde!(ContinuousScalar, X1_)
z3 = Normal(9.5,0.75)
measX1 = manikde!(ContinuousScalar, [rand(z3,1) for _ in 1:100])

# do actual product
posterioriX1 = predX1 * measX1
X1 = getPoints(posterioriX1)


## predict, 1->2 seconds

z4 = Normal(8,2.0) # odo
statemodel = LinearRelative( z4 )
X2__ = forwardConvolve(X1, statemodel)

@cast X2_[i,j] := X2__[j][i]

@test size(X2_) == (1,100)
@test 15 < Statistics.mean(X2_) < 25

##

end


@error "Add a test to ensure that approxConv does NOT change target variable VND.val values!"