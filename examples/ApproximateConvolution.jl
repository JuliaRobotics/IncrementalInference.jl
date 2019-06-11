

using IncrementalInference, KernelDensityEstimate, Distributions
using Gadfly # for draw PDF
using Test

import IncrementalInference: getSample

# switch off y ticks
toggleYTicks()


mutable struct MultiModalConditional <: IncrementalInference.FunctorPairwise
  x::Vector{Distribution}
  hypo::Categorical
  MultiModalConditional{D <: Distribution}(x::Vector{D}, p::Categorical) = new(x, p)
end
function getSample(dpl::MultiModalConditional, N::Int=1)
  d = length(dpl.hypo.p)
  p = rand(dpl.hypo, N)
  ret = zeros(1,N)
  for i in 1:N
    ret[i] = rand(dpl.x[p[i]])
  end
  return (ret, p)
end

function (dp::MultiModalConditional)(res::Vector{Float64},
                                    idx::Int,
                                    meas::Tuple{Array{Float64,2},Vector{Int64}},
                                    x1::Array{Float64},
                                    x2::Array{Float64}  )
  #
  res[1] = meas[1][1,idx] - (x2[1,idx]-x1[1,idx])
  nothing
end





fg = initfg()

N=100

doors = [-20.0, 0.0, 20.0]'
pd = kde!(doors,[2.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);
v1 = addVariable!(fg,:x1,doors,N=N)
f1  = addFactor!(fg,[v1],Obsv2( doors2, bws', [1.0])) #, samplefnc=getSample


# not initialized
v2 = addVariable!(fg,:x2, N=N)

mmc = MultiModalConditional([Normal(-5,0.5),Normal(5,0.5)],Categorical([0.5,0.5]))
f2 = addFactor!(fg, [:x1; :x2], mmc )


# Graphs.plot(fg.g)


pts = approxConv(fg, :x1x2f1, :x2)


q2 = kde!(getSample(mmc,2000)[1])
h1 = plotKDE([getVertKDE(v1), q2],c=["red";"green"],fill=true, xlbl="")
h2 = plotKDE(kde!(pts),fill=true,xlbl="", title="N = 100")



draw(PDF("approxconv.pdf",14cm,10cm),vstack(h1,h2))
# @async run(`evince approxconv.pdf`)

#
