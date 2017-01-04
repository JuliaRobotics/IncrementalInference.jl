# sqrt example

using Distributions
using KernelDensityEstimate
using IncrementalInference

import IncrementalInference: getSample



type ProductNumbers <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
getSample(s::ProductNumbers, N::Int=1) = (rand(s.z,N)', )

function (s::ProductNumbers)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      Y::Array{Float64,2},
      XY::Array{Float64,2}  )
  #
  res[1] = XY[1,idx] - X[1,idx]*Y[1,idx] + meas[1][1,idx]
  nothing
end



type Square <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
getSample(s::Square, N::Int=1) = (rand(s.z,N)', )

function (s::Square)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      XX::Array{Float64,2}  )
  #
  res[1] = XX[1,idx] - X[1,idx]*X[1,idx] + meas[1][1,idx]
  nothing
end



type AreEqual <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
function getSample(s::AreEqual, N::Int=1)
  return (rand(s.z,N)', )
end

function (s::AreEqual)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      Y::Array{Float64,2}  )
  #
  res[1] = X[1,idx]-Y[1,idx] + meas[1][1,idx]
  nothing
end


type NumbersPrior <: IncrementalInference.FunctorSingleton
  z::BallTreeDensity
end
getSample(s::NumbersPrior, N::Int=1) = (KernelDensityEstimate.sample(s.z,N)[1], )



N=300
fg = emptyFactorGraph()

x0 = [1.0+0.1*randn(100);10+0.1*randn(100)]'
addNode!(fg, :x, x0,N=N)
# addNode!(fg, :y, x0,N=N)

pts = [rand(Distributions.Normal(4.0,0.05),100);rand(Distributions.Normal(144.0,0.05),100)]
md = kde!(pts)
npx = NumbersPrior(md)
pts0 = getSample(npx,N)[1]
addNode!(fg, :xy, pts0, N=N)

addFactor!(fg, [getVert(fg, :xy)], npx)
#
# xey = AreEqual(Distributions.Normal(0.0,0.01))
# addFactor!(fg, [getVert(fg, :x);getVert(fg, :y)], xey)

# xty = ProductNumbers(Distributions.Normal(0.0,0.01))
# addFactor!(fg, [getVert(fg, :x);getVert(fg, :y);getVert(fg, :xy)], xty)


xty = Square(Distributions.Normal(0.0,0.01))
addFactor!(fg, [getVert(fg, :x);getVert(fg, :xy)], xty)



tree = wipeBuildNewTree!(fg)


@time inferOverTreeR!(fg,tree, N=N)

plotKDE(getVertKDE(fg,:xy),N=2000)

plotKDE(getVertKDE(fg,:x),N=2000)

# plotKDE(getVertKDE(fg,:y))

spyCliqMat(tree.cliques[1])


plotKDE(md)

#
