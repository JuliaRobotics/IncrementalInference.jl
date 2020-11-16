function fastnorm(u)
  # dest[1] = ...
  n = length(u)
  T = eltype(u)
  s = zero(T)
  @fastmath @inbounds @simd for i in 1:n
      s += u[i]^2
  end
  @fastmath @inbounds return sqrt(s)
end

"""
    $SIGNATURES

Sample the factor stochastic model `N::Int` times and store the samples in the preallocated `ccw.measurement` container.

DevNotes
- Use in place operations where possible and remember `measurement` is a `::Tuple`.
"""
function freshSamples(usrfnc::T, N::Int, fmd::FactorMetadata, vnd::Vector=[]) where { T <: FunctorInferenceType }
  if !hasfield(T, :specialSampler)
    getSample(usrfnc, N)
  else
    usrfnc.specialSampler(usrfnc, N, fmd, vnd...)
  end
end

function freshSamples(usrfnc::T, N::Int=1) where {T<:FunctorInferenceType}
  if hasfield(T, :specialSampler)
    error("specialSampler requires FactorMetadata and VariableNodeDatas")
  end
  freshSamples(usrfnc, N, FactorMetadata(),)
end

function freshSamples(dfg::AbstractDFG, sym::Symbol, N::Int=1)
  fct = getFactor(dfg, sym)
  usrfnc = getFactorType(fct)
  if hasfield(typeof(usrfnc), :specialSampler)
    variables = getVariable.(dfg, getVariableOrder(fct))
    fmd = _defaultFactorMetadata(variables)
    freshSamples(usrfnc, N, fmd, variables )
  else
    freshSamples(usrfnc, N)
  end
end

# TODO, add Xi::Vector{DFGVariable} if possible
function freshSamples!( ccwl::CommonConvWrapper, 
                        N::Int, 
                        fmd::FactorMetadata, 
                        vnd::Vector=[] )
  #
  # if size(ccwl.measurement, 2) == N
  # DOESNT WORK DUE TO TUPLE, not so quick and easy
  #   ccwl.measurement .= getSample(ccwl.usrfnc!, N)
  # else
    ccwl.measurement = freshSamples(ccwl.usrfnc!, N, fmd, vnd)
  # end
  nothing
end
function freshSamples!(ccwl::CommonConvWrapper, N::Int=1)
  # could maybe use default to reduce member functions
  freshSamples!(ccwl, N, FactorMetadata(),)
end

function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
  # n = length(X)
  Y = deepcopy(Alt)
  for i in 1:d
    Y[p[i]] = X[i]
  end
  return Y
end

"""
    $(SIGNATURES)

Shuffle incoming X into random positions in fr.Y.
Shuffled fr.Y will be placed back into fr.X[:,fr.gwp.particleidx] upon fr.gwp.usrfnc(x, res).
"""
function shuffleXAltD!(ccwl::CommonConvWrapper, X::Vector{Float64})
  # populate defaults from existing values
  for i in 1:ccwl.xDim
    ccwl.cpt[Threads.threadid()].Y[i] = ccwl.cpt[Threads.threadid()].X[i, ccwl.cpt[Threads.threadid()].particleidx]
  end
  # populate as many measurment dimensions randomly for calculation
  for i in 1:ccwl.zDim
    ccwl.cpt[Threads.threadid()].Y[ccwl.cpt[Threads.threadid()].p[i]] = X[i]
  end
  nothing
end


function (ccw::CommonConvWrapper)(res::AbstractVector{<:Real}, x::AbstractVector{<:Real})
  shuffleXAltD!(ccw, x)
  ccw.params[ccw.varidx][:, ccw.cpt[Threads.threadid()].particleidx] = ccw.cpt[Threads.threadid()].Y
  # evaulate the user provided residual function with constructed set of parameters
  ret = ccw.usrfnc!(res,
                    ccw.cpt[Threads.threadid()].factormetadata,
                    ccw.cpt[Threads.threadid()].particleidx,
                    ccw.measurement,
                    ccw.params[ccw.cpt[Threads.threadid()].activehypo]...) # optmize the view here, re-use the same memory
  return ret
end


function (ccw::CommonConvWrapper)(x::Vector{Float64})
  # set internal memory to that of external caller value `x`, special care if partial
  if !ccw.partial
    ccw.params[ccw.varidx][:, ccw.cpt[Threads.threadid()].particleidx] .= x #ccw.Y
  else
    ccw.params[ccw.varidx][ccw.cpt[Threads.threadid()].p, ccw.cpt[Threads.threadid()].particleidx] .= x #ccw.Y
  end
  # evaluate the user provided residual function with constructed set of parameters
  ccw.usrfnc!(ccw.cpt[Threads.threadid()].res,
              ccw.cpt[Threads.threadid()].factormetadata,
              ccw.cpt[Threads.threadid()].particleidx,
              ccw.measurement,
              ccw.params[ccw.cpt[Threads.threadid()].activehypo]...)
end



"""
    $(SIGNATURES)

Update cliq `cliqID` in Bayes (Juction) tree `bt` according to contents of `urt` -- intended use is to update main clique after a upward belief propagation computation has been completed per clique.
"""
function updateFGBT!( fg::AbstractDFG,
                      cliq::TreeClique,
                      IDvals::Dict{Symbol, TreeBelief};
                      dbg::Bool=false,
                      fillcolor::String="",
                      logger=ConsoleLogger()  )
  #
  # if dbg
  #   # TODO find better location for the debug information (this is old code)
  #   cliq.attributes["debug"] = deepcopy(urt.dbgUp)
  # end
  if fillcolor != ""
    setCliqueDrawColor!(cliq, fillcolor)
  end
  for (id,dat) in IDvals
    with_logger(logger) do
      @info "updateFGBT! up -- update $id, inferdim=$(dat.inferdim)"
    end
    updvert = DFG.getVariable(fg, id)
    setValKDE!(updvert, deepcopy(dat), true) ## TODO -- not sure if deepcopy is required
  end
  with_logger(logger) do
    @info "updateFGBT! up -- updated $(getLabel(cliq))"
  end
  nothing
end

#
