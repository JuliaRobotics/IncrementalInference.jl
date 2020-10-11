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

function numericRoot(residFnc::Function, measurement, parameters, x0::Vector{Float64})
  # function is being deprecated
  return (nlsolve(   (res, X) -> residFnc(res, measurement, parameters, X), x0, inplace=true )).zero
end

"""
    $SIGNATURES

Sample the factor stochastic model `N::Int` times and store the samples in the preallocated `ccw.measurement` container.

DevNotes
- Use in place operations where possible and remember `measurement` is a `::Tuple`.
"""
function freshSamples(usrfnc::T, N::Int, fmd::FactorMetadata, vnd::Vector=[]) where { T <: FunctorInferenceType }
  @show T, fieldnames(T)
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
    freshSamples(usrfnc, N, FactorMetadata(), getVariable.(dfg,getVariableOrder(fct)) )
  else
    freshSamples(usrfnc, N)
  end
end

# TODO, add Xi::Vector{DFGVariable} if possible
function freshSamples!(ccwl::CommonConvWrapper, N::Int, fmd::FactorMetadata, vnd::Vector=[])
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
    $SIGNATURES

Calculate both measured and predicted relative variable values, starting with `from` at zeros up to `to::Symbol`.

Notes
- assume single variable separators only.
"""
function accumulateFactorChain( dfg::AbstractDFG,
                                from::Symbol,
                                to::Symbol,
                                fsyms::Vector{Symbol}=findFactorsBetweenNaive(dfg, from, to);
                                initval=zeros(size(getVal(dfg, from))))

  # get associated variables
  svars = union(ls.(dfg, fsyms)...)

  # use subgraph copys to do calculations
  tfg_meas = buildSubgraph(dfg, [svars;fsyms])
  tfg_pred = buildSubgraph(dfg, [svars;fsyms])

  # drive variable values manually to ensure no additional stochastics are introduced.
  nextvar = from
  initManual!(tfg_meas, nextvar, initval)
  initManual!(tfg_pred, nextvar, initval)

  # nextfct = fsyms[1] # for debugging
  for nextfct in fsyms
    nextvars = setdiff(ls(tfg_meas,nextfct),[nextvar])
    @assert length(nextvars) == 1 "accumulateFactorChain requires each factor pair to separated by a single variable"
    nextvar = nextvars[1]
    meas, pred = solveFactorMeasurements(dfg, nextfct)
    pts_meas = approxConv(tfg_meas, nextfct, nextvar, (meas,ones(Int,100),collect(1:100)))
    pts_pred = approxConv(tfg_pred, nextfct, nextvar, (pred,ones(Int,100),collect(1:100)))
    initManual!(tfg_meas, nextvar, pts_meas)
    initManual!(tfg_pred, nextvar, pts_pred)
  end
  return getVal(tfg_meas,nextvar), getVal(tfg_pred,nextvar)
end


"""
    $SIGNATURES

Set the color of a cliq in the Bayes (Junction) tree.
"""
function setCliqDrawColor(cliq::TreeClique, fillcolor::String)::Nothing
  cliq.attributes["fillcolor"] = fillcolor
  cliq.attributes["style"] = "filled"
  nothing
end

function getCliqueDrawColor(cliq::TreeClique)
  haskey(cliq.attributes, "fillcolor") ? cliq.attributes["fillcolor"] : nothing
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
    setCliqDrawColor(cliq, fillcolor)
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
