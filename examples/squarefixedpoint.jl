# sqrt example

using Distributions
using KernelDensityEstimate
using IncrementalInference
using Gadfly


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


type NumbersPrior <: IncrementalInference.FunctorSingleton
  z::BallTreeDensity
end
getSample(s::NumbersPrior, N::Int=1) = (KernelDensityEstimate.sample(s.z,N)[1], )




function approxHilbertInnerProd(p::BallTreeDensity, phi::Function; N::Int=1000)
  ran = getKDERange(p)
  intrv = (ran[2]-ran[1])/N
  vec = linspace(ran..., N)
  yval = phi.(vec)
  intrv*sum(yval)
end




N=100
fg = emptyFactorGraph()
p = Dict{Symbol, BallTreeDensity}()


x0 = 0.5-rand(1,100)  #[1.0+0.1*randn(100);10+0.1*randn(100)]'
addNode!(fg, :x, x0, N=N)
# addNode!(fg, :y, x0,N=N)

pts = rand(Distributions.Normal(4.0,0.05),100)   #;rand(Distributions.Normal(144.0,0.05),100)]
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
addFactor!(fg, [:x,:xy], xty)

# Graphs.plot(fg.g)

tree = wipeBuildNewTree!(fg)
# spyCliqMat(tree.cliques[1])

iters = 100

XXh = zeros(iters, 2)
XYh = zeros(iters, 2)
XX = zeros(iters, 2)
XY = zeros(iters, 2)
XX2 = zeros(iters, 2)
XY2 = zeros(iters, 2)

coshalf = (x) -> cos(0.5*x)
sinhalf = (x) -> sin(0.5*x)
cos2 = (x) -> cos(2*x)
sin2 = (x) -> sin(2*x)


for i in 1:iters
  inferOverTree!(fg,tree, N=N)

  p[:x] = getVertKDE(fg, :x)
  p[:xy] = getVertKDE(fg, :xy)

  XXh[i,1] = approxHilbertInnerProd(p[:x], coshalf)
  XXh[i,2] = approxHilbertInnerProd(p[:x], sinhalf)
  XYh[i,1] = approxHilbertInnerProd(p[:xy], coshalf)
  XYh[i,2] = approxHilbertInnerProd(p[:xy], sinhalf)

  XX[i,1] = approxHilbertInnerProd(p[:x], cos)
  XX[i,2] = approxHilbertInnerProd(p[:x], sin)
  XY[i,1] = approxHilbertInnerProd(p[:xy], cos)
  XY[i,2] = approxHilbertInnerProd(p[:xy], sin)

  XX2[i,1] = approxHilbertInnerProd(p[:x], cos2)
  XX2[i,2] = approxHilbertInnerProd(p[:x], sin2)
  XY2[i,1] = approxHilbertInnerProd(p[:xy], cos2)
  XY2[i,2] = approxHilbertInnerProd(p[:xy], sin2)

  # plotKDE(p[:x],N=1000)
  # plotKDE(p[:xy],N=1000)
end


Gadfly.plot(x=XXh[:,1],y=XX[:,2], Geom.path())
Gadfly.plot(x=XYh[:,1],y=XY[:,2], Geom.path())

Gadfly.plot(x=XX[:,1],y=XX[:,2], Geom.path())
Gadfly.plot(x=XY[:,1],y=XY[:,2], Geom.path())

Gadfly.plot(x=XX2[:,1],y=XX[:,2], Geom.path())
Gadfly.plot(x=XY2[:,1],y=XY[:,2], Geom.path())


plotKDE(p[:x],N=1000)
plotKDE(p[:xy],N=1000)

# plotKDE(getVertKDE(fg,:y))



# Now evaluate the value function for studing the Bellman equation























#
