## Include needed packages
using IncrementalInference
using RoMEPlotting

# import getSample to be extended for user factor MultiModalConditional 
import IncrementalInference: getSample



## create a new facor type MultiModalConditional
## FIXME, this approach is unnecessary, see `::Mixture` instead.  See Caesar.jl documentation for details.
mutable struct MultiModalConditional <: AbstractRelativeRoots
  x::Vector{Distribution}
  hypo::Categorical
  MultiModalConditional(x::Vector{<:Distribution}, p::Categorical) = new(x, p)
end
function getSample(cf::CalcFactor{<:MultiModalConditional}, N::Int=1)
  d = length(cf.factor.hypo.p)
  p = rand(cf.factor.hypo, N)
  ret = zeros(1,N)
  for i in 1:N
    ret[i] = rand(cf.factor.x[p[i]])
  end
  return (ret, p)
end

function (cf::CalcFactor{<:MultiModalConditional})( res::AbstractVector{<:Real},
                                                    meas,
                                                    x1,
                                                    x2  )
  #
  res[1] = meas[1] - (x2[1]-x1[1])
  nothing
end


## build factor graph and populate
fg = initfg()

N=100

doors = [-20.0 0.0 20.0]
pd = kde!(doors,[2.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);
v1 = addVariable!(fg,:x1,ContinuousScalar,N=N)
f1 = addFactor!(fg,[v1],Prior(pd)) #, samplefnc=getSample

# not initialized
v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)

mmc = MultiModalConditional([Normal(-5,0.5),Normal(5,0.5)],Categorical([0.5,0.5]))
f2 = addFactor!(fg, [:x1; :x2], mmc )

pts = approxConv(fg, :x1x2f1, :x2)

## do some plotting
meas = sampleFactor(f2,2000)
q2 = kde!(meas)
h1 = plotKDE([getBelief(v1), q2],c=["red";"green"],fill=true, xlbl="")
h2 = plotKDE(kde!(pts),fill=true,xlbl="", title="N = 100")


draw(PDF("approxconv.pdf",14cm,10cm),vstack(h1,h2))
# @async run(`evince approxconv.pdf`)
