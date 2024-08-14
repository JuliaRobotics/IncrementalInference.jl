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
  varType::Union{InstanceType{<:InferenceVariable}, InstanceType{<:AbstractFactor}},
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
  nodeType::Union{InstanceType{<:InferenceVariable}, InstanceType{<:AbstractFactor}},
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
  _allowThreads::Bool=true
)
  #
  
  # FIXME get allocations here down to 0
    # TODO make this an in-place operation as far possible
  # TODO make this a multithreaded sampling function
  # build a CalcFactor object and get fresh samples.
  # cf = CalcFactor(ccwl; _allowThreads) 
  resize!(ccwl.measurement, N)
  ccwl.measurement[:] = sampleFactor(ccwl, N; _allowThreads)

  return ccwl.measurement
end

function sampleFactor(
  ccwl::CommonConvWrapper, 
  N::Int; 
  _allowThreads::Bool=true
)
  #
  cf = CalcFactorNormSq(ccwl; _allowThreads) 
  return sampleFactor(cf, N)
end

sampleFactor(
  fct::DFGFactor, 
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
  fct::AbstractFactor,
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

"""
    $SIGNATURES

Check if a variable might already be located at the test location, by means of a (default) `refKey=:simulated` PPE stored in the existing variables.

Notes
- Checks, using provided `factor` from `srcLabel` in `fg` to an assumed `dest` variable whcih may or may not yet exist.
- This function was written to aid in building simulation code, 
  - it's use in real world usage may have unexpected behaviour -- hence not exported.
- Return `::Tuple{Bool, Vector{Float64}, Symbol}`, eg. already exists `(true, [refVal], :l17)`, or if a refernce variable does not yet `(false, [refVal], :l28)`.
  - Vector contains the PPE reference location of the new variable as calculated.
- Auto `destPrefix` is trying to parse `destRegex` labels like `l\\d+` or `tag\\d+`, won't work with weirder labels e.g. `:l_4_23`.
  - User can overcome weird names by self defining `destPrefix` and `srcNumber`.
  - User can also ignore and replace the generated new label `Symbol(destPrefix, srcNumber)`.
- This function does not add new variables or factors to `fg`, user must do that themselves after.
  - Useful to use in combination with `setPPE!` on new variable.
- At time of writing `accumulateFactorMeans` could only incorporate priors or binary relative factors.
  - internal info, see [`solveFactorParametric`](@ref),
  - This means at time of writing `factor` must be a binary factor.
- Tip, if simulations are inducing odometry bias, think of using two factors from caller (e.g. simPerfect and simBias).

Example
```julia
# fg has :x5 and :l2 and PPEs :simulated exists in all variables
# user wants to add a factor from :x5 to potential new :l5, but maybe a (simulated) variable, say :l2, is already there.

newFactor = RoME.Pose2Point2BearingRange(Normal(), Normal(20,0.5))
isAlready, simPPE, genLabel = IIF._checkVariableByReference(fg, :x5, r"l\\d+", Point2, newFactor)

# maybe add new variable
if !isAlready
  @info "New variable with simPPE" genLabel simPPE 
  newVar = addVariable!(fg, genLabel, Point2)
  addFactor!(fg, [:x5; genLabel], newFactor)

  # also set :simulated PPE for similar future usage
  newPPE = DFG.MeanMaxPPE(:simulated, simPPE, simPPE, simPPE)
  setPPE!(newVar, :simulated, typeof(newPPE), newPPE)   # TODO this API can be improved
else
  @info "Adding simulated loop closure with perfect data association" :x5 genLabel
  addFactor!(fg, [:x5; genLabel], newFactor)
end

# the point is that only the (0,20) values in newFactor are needed, all calculations are abstracted away.
```

See also: [`RoME.generateGraph_Honeycomb!`](@ref), [`accumulateFactorMeans`](@ref), [`getPPE`](@ref)
"""
function _checkVariableByReference(
  fg::AbstractDFG,
  srcLabel::Symbol,
  destRegex::Regex,
  destType::Type{<:InferenceVariable},
  factor::AbstractRelative;
  srcType::Type{<:InferenceVariable} = getVariableType(fg, srcLabel) |> typeof,
  doRef::Bool = true,
  refKey::Symbol = :simulated,
  prior = if !doRef
    nothing
  else
    DFG._getPriorType(srcType)(
    MvNormal(getPPE(fg[srcLabel], refKey).suggested, diagm(ones(getDimension(srcType)))),
  )
  end,
  atol::Real = 1e-2,
  destPrefix::Symbol = match(r"[a-zA-Z_]+", destRegex.pattern).match |> Symbol,
  srcNumber = match(r"\d+", string(srcLabel)).match |> x -> parse(Int, x),
  overridePPE = doRef ? nothing : zeros(getDimension(destType)),
)
  #

  refVal = if overridePPE !== nothing
    overridePPE
  else
    # calculate and add the reference value
    # TODO refactor consolidation to use `_buildGraphByFactorAndTypes!`
    tfg = initfg()
    addVariable!(tfg, :x0, srcType)
    addFactor!(tfg, [:x0], prior)
    addVariable!(tfg, :l0, destType)
    addFactor!(tfg, [:x0; :l0], factor; graphinit = false)

    # calculate where the landmark reference position is
    accumulateFactorMeans(tfg, [:x0f1; :x0l0f1])
  end

  ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)

  # now check if we already have a landmark at this location
  varLms = ls(fg, destRegex) |> sortDFG
  already = if doRef
    ppeLms = getPPE.(getVariable.(fg, varLms), refKey) .|> x -> x.suggested
    errmask = ppeLms .|> (x -> isapprox(x, refVal; atol = atol))
    any(errmask)
  else
    false
  end

  if already
    # does exist, ppe, variableLabel
    alrLm = varLms[findfirst(errmask)]
    # @info "Variable on :$refKey does exists at" srcLabel alrLm
    return true, ppe, alrLm
  end

  # Nope does not exist, ppe, generated new variable label only
  return false, ppe, Symbol(destPrefix, srcNumber)
end

function _checkVariableByReference(
  fg::AbstractDFG,
  srcLabel::Symbol,
  destRegex::Regex,
  destType::Type{<:InferenceVariable},
  factor::AbstractPrior;
  srcType::Type{<:InferenceVariable} = getVariableType(fg, srcLabel) |> typeof,
  doRef::Bool = true,
  refKey::Symbol = :simulated,
  prior = typeof(factor)(MvNormal(getMeasurementParametric(factor)...)),
  atol::Real = 1e-3,
  destPrefix::Symbol = match(r"[a-zA-Z_]+", destRegex.pattern).match |> Symbol,
  srcNumber = match(r"\d+", string(srcLabel)).match |> x -> parse(Int, x),
  overridePPE = doRef ? nothing : zeros(getDimension(destType)),
)
  #

  refVal = if overridePPE !== nothing
    overridePPE
  else
    getMeasurementParametric(factor)[1]
  end

  ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)

  # Nope does not exist, ppe, generated new variable label only
  return false, ppe, Symbol(destPrefix, srcNumber)
end

#
