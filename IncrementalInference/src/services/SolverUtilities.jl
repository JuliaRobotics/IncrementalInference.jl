function fastnorm(u)
  # dest[1] = ...
  n = length(u)
  T = eltype(u)
  s = zero(T)
  @fastmath @inbounds @simd for i = 1:n
    s += u[i]^2
  end
  @fastmath @inbounds return sqrt(s)
end

# """
#     $TYPEDSIGNATURES

# Calculate the Kernel Embedding MMD 'distance' between sample points (or kernel density estimates).

# Notes
# - `bw::Vector=[0.001;]` controls the mmd kernel bandwidths.
# - Overloading from ApproxManifoldProducts

# Related

# `AMP.kld`
# """
function mmd(
  p1::AbstractVector{P1},
  p2::AbstractVector{P2},
  varType::Union{InstanceType{<:StateType}, InstanceType{<:AbstractObservation}},
  threads::Bool = true;
  bw::AbstractVector{<:Real} = SA[0.001;],
) where {P1 <: AbstractVector, P2 <: AbstractVector}
  #
  mani = getManifold(varType)
  return mmd(mani, p1, p2, length(p1), length(p2), threads; bw)
end

function mmd(
  p1::ManifoldKernelDensity,
  p2::ManifoldKernelDensity,
  nodeType::Union{InstanceType{<:StateType}, InstanceType{<:AbstractObservation}},
  threads::Bool = true;
  bw::AbstractVector{<:Real} = SA[0.001;],
  asPartial::Bool = true
)
  #
  return mmd(getPoints(p1, asPartial), getPoints(p2, asPartial), nodeType, threads; bw)
end

# part of consolidation, see #927
function sampleFactor!(
  ccwl::CommonConvWrapper, 
  N::Int; 
  _allowThreads::Bool=true,
  keepCalcFactor::Union{Nothing, <:Channel} = nothing,
)
  #
  
  # FIXME get allocations here down to 0
    # TODO make this an in-place operation as far possible
  # TODO make this a multithreaded sampling function
  # build a CalcFactor object and get fresh samples.
  # cf = CalcFactor(ccwl; _allowThreads) 
  resize!(ccwl.measurement, N)
  ccwl.measurement[:] = sampleFactor(ccwl, N; _allowThreads, keepCalcFactor)

  return ccwl.measurement
end

function sampleFactor(
  ccwl::CommonConvWrapper, 
  N::Int; 
  _allowThreads::Bool=true,
  keepCalcFactor::Union{Nothing, <:Channel} = nothing,
)
  #
  cf = CalcFactorNormSq(ccwl; _allowThreads) 
  smpls = sampleFactor(cf, N)
  isnothing(keepCalcFactor) ? nothing : put!(keepCalcFactor, cf)
  return smpls 
end

sampleFactor(
  fct::FactorCompute, 
  N::Int = 1; 
  _allowThreads::Bool=true
) = sampleFactor(
  _getCCW(fct), 
  N; 
  _allowThreads
)

function sampleFactor(
  dfg::AbstractDFG, 
  sym::Symbol, 
  N::Int = 1; 
  _allowThreads::Bool=true
)
  #
  return sampleFactor(getFactor(dfg, sym), N; _allowThreads)
end

"""
    $(SIGNATURES)

Update cliq `cliqID` in Bayes (Juction) tree `bt` according to contents of `urt`.
Intended use is to update main clique after a upward belief propagation computation 
has been completed per clique.
"""
function updateFGBT!(
  fg::AbstractDFG,
  cliq::TreeClique,
  IDvals::Dict{Symbol, TreeBelief};
  dbg::Bool = false,
  fillcolor::String = "",
  logger = ConsoleLogger(),
)
  #
  # if dbg
  #   # TODO find better location for the debug information (this is old code)
  #   cliq.attributes["debug"] = deepcopy(urt.dbgUp)
  # end
  if fillcolor != ""
    setCliqueDrawColor!(cliq, fillcolor)
  end
  for (id, dat) in IDvals
    with_logger(logger) do
      @info "updateFGBT! up -- update $id, infoPerCoord=$(dat.infoPerCoord)"
    end
    updvert = DFG.getVariable(fg, id)
    setValKDE!(updvert, deepcopy(dat), true) ## TODO -- not sure if deepcopy is required
  end
  with_logger(logger) do
    @info "updateFGBT! up -- updated $(getLabel(cliq))"
  end
  return nothing
end

"""
    $SIGNATURES

Build a graph given one factor and an ordered vector of `(variables types,nothing).  In addition, init values can be passed instead of nothing.

Notes
- Often used to quickly generate temporary graphs for a variety of local calculations.
- does not yet support split `_` characters in auto-find `lastVar` from `varPattern`. 
- Will always add a factor, but will skip adding variable labels that already exist in `dfg`.

DevNotes
- TODO allow pts to be full MKD beliefs, part of replacing old `approxConvCircular`, see #1351
"""
function _buildGraphByFactorAndTypes!(
  fct::AbstractObservation,
  varTypes::Tuple,
  pts::Tuple = ();
  dfg::AbstractDFG = initfg(),
  solveKey::Symbol = :default,
  newFactor::Bool = true,
  destPattern::Regex = r"x\d+",
  destPrefix::Symbol = match(r"[a-zA-Z_]+", destPattern.pattern).match |> Symbol,
  _allVars::AbstractVector{Symbol} = sortDFG(ls(dfg, destPattern)),
  currLabel::Symbol = 0 < length(_allVars) ? _allVars[end] : Symbol(destPrefix, 0),
  currNumber::Integer = reverse(match(r"\d+", reverse(string(currLabel))).match) |>
                        x -> parse(Int, x),
  graphinit::Bool = false,
  _blockRecursion::Bool = false,
)
  #

  # TODO generalize beyond binary
  len = length(varTypes)
  vars = Symbol[Symbol(destPrefix, s_) for s_ in (currNumber .+ (1:len))]
  for (s_, vTyp) in enumerate(varTypes)
    # add the necessary variables
    exists(dfg, vars[s_]) ? nothing : addVariable!(dfg, vars[s_], vTyp)
    # set the numerical values if available
    # TODO allow pts to come in as full MKD beliefs, not just one point
    if ((0 < length(pts)) && (pts[s_] isa Nothing))
      nothing
    else
      initVariable!(dfg, vars[s_], [pts[s_]], solveKey; bw = ones(getDimension(vTyp)))
    end
  end
  # if newFactor then add the factor on vars, else assume only one existing factor between vars
  _dfgfct = if newFactor
    addFactor!(dfg, vars, fct; graphinit, _blockRecursion)
  else
    getFactor(dfg, intersect((ls.(dfg, vars))...)[1])
  end

  return dfg, _dfgfct
end

#
