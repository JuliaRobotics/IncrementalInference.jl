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
    $TYPEDSIGNATURES

Calculate the Kernel Embedding MMD 'distance' between sample points (or kernel density estimates).

Notes
- `bw::Vector=[0.001;]` controls the mmd kernel bandwidths.

Related

`KDE.kld`
"""
function mmd( p1::AbstractMatrix{<:Real}, 
              p2::AbstractMatrix{<:Real}, 
              varType::Union{InstanceType{InferenceVariable},InstanceType{FunctorInferenceType}};
              bw::AbstractVector{<:Real}=[0.001;] )
  #
  manis = convert(AMP.Manifold, varType)
  mmd(p1, p2, manis, bw=bw)  
end

function mmd( p1::BallTreeDensity, 
              p2::BallTreeDensity, 
              nodeType::Union{InstanceType{InferenceVariable},InstanceType{FunctorInferenceType}};
              bw::AbstractVector{<:Real}=[0.001;])
  #
  mmd(getPoints(p1), getPoints(p2), nodeType, bw=bw)
end

# moved to CalcFactor.jl


function freshSamples(ccwl::CommonConvWrapper,
                      N::Int  )
  #
  cf = CalcFactor( ccwl.usrfnc!, _getFMdThread(ccwl), 0, length(ccwl.measurement), ccwl.measurement, ccwl.params)
  freshSamples(cf, N)
end

# part of consolidation, see #927
function freshSamples!( ccwl::CommonConvWrapper, 
                        N::Int, 
                        fmd::FactorMetadata=_getFMdThread(ccwl), 
                        vnd=nothing )
  #
  # depr warning added before IIF v0.20
  vnd !== nothing ? @warn("freshSamples! no longer accepts vnd::Vector as meaningful input.") : nothing
  
  # build a CalcFactor object and get fresh samples.
  cf = CalcFactor( ccwl.usrfnc!, fmd, 0, length(ccwl.measurement), ccwl.measurement, ccwl.params)
  # TODO make this an in-place operation as far possible
  ccwl.measurement = freshSamples(cf, N)    
  nothing
end


function freshSamples(dfg::AbstractDFG, 
                      sym::Symbol, 
                      N::Int=1 )
  #
  # fct = getFactor(dfg, sym)
  return freshSamples(_getCCW(dfg, sym), N)
  # usrfnc = getFactorType(fct)
  # variables = getVariable.(dfg, getVariableOrder(fct))
  
  # FIXME why not use ccwl.cpt[thrid].fmd rather than building a new one 
  # fmd = _getFMdThread() # what about shared mem concurrancy?
  # FIXME fix/avoid getSample issue in testMultihypoFMD.jl: ummm, can we sample without knowing the hypotheses?
  # also deconv?
  # not really, because that would imply stochastic dependency on association before noise process??
  # fmd = FactorMetadata(variables, getLabel.(variables), Vector{Matrix{Float64}}(), :null, nothing)
  # freshSamples(usrfnc, N, fmd )
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
