# sqrt example

using Distributions
using KernelDensityEstimate, KernelDensityEstimatePlotting
using IncrementalInference
using Gadfly, DataFrames


import IncrementalInference: getSample



mutable struct ProductNumbers <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
getSample(s::ProductNumbers, N::Int=1) = (rand(s.z,N)', )

function (s::ProductNumbers)(res::Array{Float64},
      userdata::FactorMetadata,
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      Y::Array{Float64,2},
      XY::Array{Float64,2}  )
  #
  res[1] = XY[1,idx] - X[1,idx]*Y[1,idx] + meas[1][1,idx]
  nothing
end



struct AreEqual <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
function getSample(s::AreEqual, N::Int=1)
  return (rand(s.z,N)', )
end

function (s::Square)(res::Array{Float64},
      userdata::FactorMetadata,
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      Y::Array{Float64,2}  )
  #
  res[1] = X[1,idx]-Y[1,idx] + meas[1][1,idx]
  nothing
end



struct Square <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
getSample(s::Square, N::Int=1) = (rand(s.z,N)', )

function (s::AreEqual)(res::Array{Float64},
      userdata::FactorMetadata,
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      XX::Array{Float64,2}  )
  #
  res[1] = XX[1,idx] - X[1,idx]*X[1,idx] + meas[1][1,idx]
  nothing
end


mutable struct NumbersPrior <: IncrementalInference.FunctorSingleton
  z::BallTreeDensity
end
getSample(s::NumbersPrior, N::Int=1) = (KernelDensityEstimate.sample(s.z,N)[1], )




function extractCycleProjections(FGl::Vector{FactorGraph}, FR; N::Int=1000)
  itersl = length(FGl)
  numfr = length(FR)
  ALLXX = zeros(numfr, itersl, 2)
  ALLXY = zeros(numfr, itersl, 2)
  PL = Array{Any,2}(itersl,3)
  p = Dict{Symbol, BallTreeDensity}()
  PP = Dict{Symbol, Vector{BallTreeDensity}}()
  DIV = Dict{Symbol, Vector{Float64}}()
  DIVREF = Dict{Symbol, Vector{Float64}}()
  PP[:x] = BallTreeDensity[]
  PP[:xy] = BallTreeDensity[]
  DIV[:x] = Float64[]
  DIV[:xy] = Float64[]
  DIVREF[:x] = Float64[]
  DIVREF[:xy] = Float64[]

  for i in 1:itersl
    p[:x] = getVertKDE(FGl[i], :x)
    p[:xy] = getVertKDE(FGl[i], :xy)

    push!(PP[:x], deepcopy(p[:x]))
    push!(PP[:xy], deepcopy(p[:xy]))

    for j in 1:numfr
      ALLXX[j,i,1] = approxHilbertInnerProd(p[:x], (x) -> phic(x, f=FR[j]), N=N)
      ALLXX[j,i,2] = approxHilbertInnerProd(p[:x], (x) -> phis(x, f=FR[j]), N=N)
      ALLXY[j,i,1] = approxHilbertInnerProd(p[:xy], (x) -> phic(x, f=FR[j]), N=N)
      ALLXY[j,i,2] = approxHilbertInnerProd(p[:xy], (x) -> phis(x, f=FR[j]), N=N)
    end

    PL[i,1] = plotKDE(p[:x],N=N, extend=0.2)
    PL[i,2] = plotKDE(p[:xy],N=N, extend=0.2)
    PL[i,3] = plotKDE([p[:x],p[:xy]],c=["red","green"],N=N, extend=0.2)

    if i > 1
      # KL-divergence
      mkl = minimum([kld(p[:x],getVertKDE(FGl[i-1],:x))[1]; kld(getVertKDE(FGl[i-1],:x),p[:x])[1]])
      push!(DIV[:x], abs(mkl))
      mklxy = minimum([kld(p[:xy],getVertKDE(FGl[i-1],:xy))[1]; kld(getVertKDE(FGl[i-1],:xy),p[:xy])[1]])
      push!(DIV[:xy], abs(mklxy))
    end

    # reference KL-divergence
    mkl = minimum([kld(p[:x],getVertKDE(FGl[end],:x))[1]; kld(getVertKDE(FGl[end],:x),p[:x])[1]])
    push!(DIVREF[:x], abs(mkl))
    mklxy = minimum([kld(p[:xy],getVertKDE(FGl[end],:xy))[1]; kld(getVertKDE(FGl[end],:xy),p[:xy])[1]])
    push!(DIVREF[:xy], abs(mklxy))
  end

  cycle = Dict{Symbol, Any}()
  # save everything
  cycle[:N] = N
  cycle[:iters] = itersl-1
  cycle[:FR] = FR
  cycle[:numfr] = numfr
  cycle[:PP] = PP
  cycle[:ALLXX] = ALLXX
  cycle[:ALLXY] = ALLXY
  cycle[:DIV] = DIV
  cycle[:DIVREF] = DIVREF
  return cycle
end



function runFullBatchIterations(;N=100, iters=50)
  fg = initfg()

  x0 = 0.5-rand(1,N)  #[1.0+0.1*randn(N);10+0.1*randn(N)]'
  addVariable!(fg, :x, x0, N=N)
  # addVariable!(fg, :y, x0,N=N)

  pts = rand(Distributions.Normal(4.0,0.05),N)   #;rand(Distributions.Normal(144.0,0.05),N)]
  md = kde!(pts)
  npx = NumbersPrior(md)
  pts0 = getSample(npx,N)[1]
  addVariable!(fg, :xy, pts0, N=N)

  addFactor!(fg, [getVariable(fg, :xy)], npx)
  #
  # xey = AreEqual(Distributions.Normal(0.0,0.01))
  # addFactor!(fg, [getVariable(fg, :x);getVariable(fg, :y)], xey)

  # xty = ProductNumbers(Distributions.Normal(0.0,0.01))
  # addFactor!(fg, [getVariable(fg, :x);getVariable(fg, :y);getVariable(fg, :xy)], xty)

  xty = Square(Distributions.Normal(0.0,0.01))
  addFactor!(fg, [:x,:xy], xty)
  # Graphs.plot(fg.g)
  tree = wipeBuildNewTree!(fg)
  # spyCliqMat(tree.cliques[1])

  FG = Vector{FactorGraph}(iters+1)

  FG[1] = deepcopy(fg)
  for i in 1:iters
    inferOverTree!(fg, tree, N=N)
    FG[i+1] = deepcopy(fg)
  end
  return FG
end


# Do all the runs
# frequencies of interest

FR = range(0.5/(2pi),stop=3.0/(2pi), length=8)
# FR = range(0.5/(2pi),stop=3.0/(2pi), length=8)


mc = 3

# data containers
FG = Vector{Vector{FactorGraph}}(undef, mc)
CYCLE = Vector{Dict}(undef, mc)


for i in 1:mc
  FG[i] = runFullBatchIterations(;N=100, iters=50)
  CYCLE[i] = extractCycleProjections(FG[i], FR, N=2000)
end





# Analyse the data

function plotXandXYFreq(CYCLE, FR, fridx)
  XXfrs = DataFrame[]
  for i in 1:mc
    push!(XXfrs, DataFrame(
      x=CYCLE[i][:ALLXX][fridx,:,1],
      y=CYCLE[i][:ALLXX][fridx,:,2],
      MC="$(i)"
    ))
  end

  XXrepeat = vcat(XXfrs...)

  plxrep = Gadfly.plot(XXrepeat,
    Geom.path(),
    Guide.title("X"),
    Guide.title("μ=0, frequency=$(round(FR[fridx],3))"),
    Guide.ylabel("p ⋅ ϕs"),
    Guide.xlabel("p ⋅ ϕc"),
    x=:x,
    y=:y,
    color=:MC
  )

  XYfrs = DataFrame[]
  for i in 1:mc
    push!(XYfrs, DataFrame(
      x=CYCLE[i][:ALLXY][fridx,:,1],
      y=CYCLE[i][:ALLXY][fridx,:,2],
      MC="$(i)"
    ))
  end

  XYrepeat = vcat(XYfrs...)

  plxyrep = Gadfly.plot(XYrepeat,
    Geom.path(),
    Guide.title("X^2"),
    Guide.title("μ=0, frequency=$(round(FR[fridx],3))"),
    Guide.ylabel("p ⋅ ϕs"),
    Guide.xlabel("p ⋅ ϕc"),
    x=:x,
    y=:y,
    color=:MC
  )
  return plxyrep,plxrep
end


for i in 1:length(FR)
  plxyrep,plxrep = plotXandXYFreq(CYCLE, FR, i)
  Gadfly.draw(PDF("sqrtexamplexxtrajFR$(i).pdf",15cm,7cm),hstack(plxyrep,plxrep))
end
# @async run(`evince sqrtexamplexxtrajFR2.pdf`)



0


#
#
# plotKDE(
# [CYCLE[1][:PP][:x][1];
# CYCLE[1][:PP][:x][2];
# CYCLE[1][:PP][:x][3]],
# c=["black","red","green"]
# )





DFs = DataFrame[]

for i in [1,2,3,5,13]
  p = CYCLE[1][:PP][:x][i]
  mxmx = getKDERange(p)
  x = [range(mxmx[1], stop=mxmx[2], length=2000);]
  push!(DFs, DataFrame(
    x = x,
    y = clamp(evaluateDualTree(p,x),0,4),
    Iteration="$(i-1)"
  ))
end

plx = Gadfly.plot(vcat(DFs...) , x=:x, y=:y, color=:Iteration,
  Geom.line,
  Guide.xlabel("X"),
  Guide.ylabel("pdf")
)


DFs = DataFrame[]

for i in [1,2,3,5,13]
  p = CYCLE[1][:PP][:xy][i]
  mxmx = getKDERange(p, extend=0.4)
  x = [range(mxmx[1], stop=mxmx[2], length=2000);]
  push!(DFs, DataFrame(
    x = x,
    y = clamp(evaluateDualTree(p,x),0,6),
    Iteration="$(i-1)"
  ))
end


plxy = Gadfly.plot(vcat(DFs...) , x=:x, y=:y, color=:Iteration,
  Geom.line,
  Guide.xlabel("X^2"),
  Guide.ylabel("pdf")
)


plh = vstack(plxy, plx)

Gadfly.draw(PDF("sqrtexamplebeliefs.pdf",12cm,9cm),plh)

@async run(`evince sqrtexamplebeliefs.pdf`)


0



# Also plot the KL divergences



DFs = DataFrame[]

for i in 1:mc
  push!(DFs, DataFrame(
    x = 1:length(CYCLE[i][:DIV][:xy]),
    y = CYCLE[i][:DIV][:xy],
    MC="$(i), X^2"
  ))
end

for i in 1:mc
  push!(DFs, DataFrame(
    x = 1:length(CYCLE[i][:DIV][:x]),
    y = CYCLE[i][:DIV][:x],
    MC="$(i), X"
  ))
end

pldiv = Gadfly.plot(vcat(DFs...) , x=:x, y=:y, color=:MC,
  Geom.line,
  Guide.title("Relative between iterations"),
  Guide.xlabel("Iterations"),
  Guide.ylabel("KL-Divergence")
)


Gadfly.draw(PDF("sqrtexamplekldrelative.pdf",10cm,6cm),pldiv)
# @async run(`evince sqrtexamplekldrelative.pdf`)




DFs = DataFrame[]
# x0 = [1.0+0.1*randn(100);10+0.1*randn(100)]'
# addVariable!(fg, :x, ContinuousScalar, N=N)
#
# pts = [rand(Distributions.Normal(4.0,0.05),100);rand(Distributions.Normal(144.0,0.05),100)]
# md = kde!(pts)
# npx = NumbersPrior(md)
# pts0 = getSample(npx,N)[1]
# addVariable!(fg, :xy, ContinuousScalar, N=N)


for i in 1:mc
  push!(DFs, DataFrame(
    x = 1:length(CYCLE[i][:DIVREF][:xy]),
    y = CYCLE[i][:DIVREF][:xy],
    MC="$(i), X^2"
  ))
end

for i in 1:mc
  push!(DFs, DataFrame(
    x = 1:length(CYCLE[i][:DIVREF][:x]),
    y = CYCLE[i][:DIVREF][:x],
    MC="$(i), X"
  ))
end

pldivref = Gadfly.plot(vcat(DFs...) , x=:x, y=:y, color=:MC,
  Geom.line,
  Guide.title("Referenced to final belief"),
  Guide.xlabel("Iterations"),
  Guide.ylabel("KL-Divergence")
)

Gadfly.draw(PDF("sqrtexamplekldreferenced.pdf",10cm,6cm),pldivref)
@async run(`evince sqrtexamplekldreferenced.pdf`)




0


#
# Gadfly.plot(x=cycle1[:ALLXX][1,:,1],y=cycle1[:ALLXX][1,:,2],
# Geom.path(), Guide.title("μ=0, frequency=$(round(FR[1],3))"))
# Gadfly.plot(x=ALLXY[1,:,1],y=ALLXY[1,:,2],
# Geom.path(), Guide.title("μ=0, frequency=$(round(FR[1],3))"))
#
# Gadfly.plot(x=ALLXX[2,:,1],y=ALLXX[2,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[2],3))"))
# Gadfly.plot(x=ALLXY[2,:,1],y=ALLXY[2,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[2],3))"))
#
# Gadfly.plot(x=ALLXX[3,:,1],y=ALLXX[2,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[3],3))"))
# Gadfly.plot(x=ALLXY[3,:,1],y=ALLXY[2,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[3],3))"))
#
# Gadfly.plot(x=ALLXX[4,:,1],y=ALLXX[2,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[4],3))"))
# Gadfly.plot(x=ALLXY[4,:,1],y=ALLXY[2,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[4],3))"))
#
# Gadfly.plot(x=ALLXX[end,:,1],y=ALLXX[3,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[5],3))"))
# Gadfly.plot(x=ALLXY[end,:,1],y=ALLXY[3,:,2], Geom.path(), Guide.title("μ=0, frequency=$(round(FR[5],3))"))
#
#
#
# # plotKDE(getVertKDE(fg,:y))
#
# PL[1,3]
# PL[2,3]
# PL[3,3]
# PL[4,3]
# PL[6,3]
#
# PL[8,3]
#
# PL[10,3]
#
# PL[12,3]
#
# PL[20,3]
#




# Now evaluate the value function for studing the Bellman equation









using KernelDensityEstimatePlotting


plot(getVertKDE(fg,:xy),N=2000)

plot(getVertKDE(fg,:x),N=2000)

# plot(getVertKDE(fg,:y))



plot(md)

#
