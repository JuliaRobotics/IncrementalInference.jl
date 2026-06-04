# Solver type structs for dispatch
struct NLLSSolver end
@kwdef struct NPBPSolver
  defaultNumKernels::Int = 128
end

"""
    $SIGNATURES

Create an uninitialized `State` on variable `v` with `num_kernels` identity-element particles.

See also: [`prepare!`](@ref), [`initAll!`](@ref)
"""
function prepareState!(
  v::VariableCompute,
  solver::Union{<:NPBPSolver, <:NLLSSolver},
  statelabel::Symbol;
  num_kernels::Int = solver.defaultNumKernels, # ensure 1 for parametric case
  varType::StateType = DFG.getStateKind(v),
)
  # work around during consolidation
  _nkrs = solver isa NLLSSolver ? 1 : num_kernels
  
  hasState(v, statelabel) && return 0
  dims = getDimension(v)
  @assert getPointType(varType) != DataType "cannot add manifold point type $(getPointType(varType)), make sure the identity element argument in @defStateType $varType arguments is correct"
  ϵ = getPointIdentity(varType)
  belief = HomotopyDensity_legacy(
    varType, 
    [ϵ for _ in 1:_nkrs]; 
    bw = zeros(dims), 
    newbw = false
  )
  mergeState!(
    v,
    State(statelabel, varType; belief, initialized = false, marginalized = false)
  )
  return 1
end
  # TODO review and refactor this function, exists as legacy from pre-v0.3.0
  # this should be the only function allocating memory for the node points (unless number of points are changed)
  # (val, bw) = if initialized
  #   pN = getBelief(v)
  #   bw = getBW(pN)[1]
  #   pNpts = getPoints(pN)
  #   isinit = true
  #   (pNpts, bw)
  # else
  #   sp = round.(Int, range(dodims; stop = dodims + dims - 1, length = dims))
  #   @assert getPointType(varType) != DataType "cannot add manifold point type $(getPointType(varType)), make sure the identity element argument in @defStateType $varType arguments is correct"
  #   val = Vector{getPointType(varType)}(undef, N)
  #   for i = 1:length(val)
  #     val[i] = getPointIdentity(varType)
  #   end
  #   bw = diagm(ones(dims))
  #   #
  #   (val, bw)
  # end

"""
    $SIGNATURES

Ensure all variables in `dfg` (matching `whereSolvable`) have a `State` for `statelabel`.
Returns the number of states newly created (0 for variables that already had the state).

See also: [`prepareState!`](@ref), [`prepare!`](@ref)
"""
function prepareStates!(
  dfg::AbstractDFG,
  solver::NPBPSolver,
  statelabel::Symbol;
  num_kernels::Int = solver.defaultNumKernels,
  whereSolvable = >=(1),
)
  count = 0
  for vl in listVariables(dfg; whereSolvable)
    vari = getVariable(dfg, vl)
    count += prepareState!(vari, solver, statelabel; num_kernels)
  end
  return count
end

function prepareStates!(
  dfg::AbstractDFG,
  solver::NLLSSolver,
  statelabel::Symbol;
  whereSolvable = >=(1),
)
  count = 0
  for vl in listVariables(dfg; whereSolvable)
    vari = getVariable(dfg, vl)
    count += prepareState!(vari, solver, statelabel)
  end
  return count
end

"""
    $SIGNATURES

For variables in `varList` check and if necessary make solverData objects for both `:default` and `:parametric` solveKeys. 

Example
```julia
num_made = makeSolverData(fg; solveKey=:parametric)
```

Notes
- Part of solving JuliaRobotics/IncrementalInference.jl issue 1637

DevNotes
- TODO, assumes parametric solves will always just be in solveKey `:parametric`.

See also: [`doautoinit!`](@ref), [`initAll!`](@ref)
"""
function makeSolverData!(
  dfg::AbstractDFG;
  solvable = 1,
  varList::AbstractVector{Symbol} = ls(dfg; whereSolvable = >=(solvable)),
  solveKey::Symbol=:default,
  N::Int = getSolverParams(dfg).N,
)
  Base.depwarn("`makeSolverData!` is deprecated, use `prepareStates!(dfg, solver, solveKey)` instead.", :makeSolverData!)
  count = 0
  for vl in varList
    v = getVariable(dfg,vl)
    if solveKey != :parametric
        count += prepareState!(v, NPBPSolver(), solveKey; num_kernels=N)
    else
        count += prepareState!(v, NLLSSolver(), solveKey)
    end
  end

  return count
end

"""
    $SIGNATURES

Return `(::Bool, ::OKVarlist, ::NotOkayVarList)` on whether all other variables (besides `loovar::Symbol`)
attached to factor `fct::Symbol` are all initialized -- i.e. `fct` is usable.

Notes:
- Special carve out for multihypo cases, see issue 427, where at least one hypothesis should be available, but not all required at first.

Development Notes
* TODO get faster version of isInitialized for database version

Related

doautoinit!, initVariable!, isInitialized, isMultihypo
"""
function factorCanInitFromOtherVars(
  dfg::AbstractDFG,
  fct::Symbol,
  loovar::Symbol;
  solveKey::Symbol = :default,
)
  #
  # all variables attached to this factor
  varsyms = listNeighbors(dfg, fct)

  # which element is being solved for
  sfidx = (1:length(varsyms))[varsyms .== loovar][1]
  # list of factors to use in init operation
  fctlist = Symbol[]
  # list fo variables that cannot be used
  faillist = Symbol[]
  isinit = Bool[]
  for vsym in varsyms
    # check each variable one by one
    xi = DFG.getVariable(dfg, vsym)
    isi = hasState(xi, solveKey) && isInitialized(xi, solveKey)
    push!(isinit, isi)
    if !isi
      push!(faillist, vsym)
    end
  end

  ## determine if this factor can be used
  # priors and general n-ary cases
  canuse = length(varsyms) == 1 || (length(faillist) == 1 && loovar in faillist)
  ## special multihypo case (at least one hypothesis is available or initializing first hypo)
  fctnode = getFactor(dfg, fct)
  # @show canuse, isMultihypo(fctnode), isinit
  if !canuse && isMultihypo(fctnode)
    # multihypo=[1;0.5;0.5] : sfidx=1, isinit=[0,1,0] -- true
    # multihypo=[1;0.5;0.5] : sfidx=1, isinit=[0,0,1] -- true
    # multihypo=[1;0.5;0.5] : sfidx=2|3, isinit=[1,0,0] -- true
    mhp = getMultihypoDistribution(fctnode).p
    allmhp, certainidx, uncertnidx = getHypothesesVectors(mhp)
    if isLeastOneHypoAvailable(sfidx, certainidx, uncertnidx, isinit)
      # special case works
      @info "allowing init from incomplete set of previously initialized hypotheses, fct=$fct"
      canuse = true
    end
  end

  # should add the factor for use?
  if canuse
    push!(fctlist, fct)
  end

  # return if can use, the factor in an array, and the non-initialized variables attached to the factor
  return (canuse, fctlist, faillist)::Tuple{Bool, Vector{Symbol}, Vector{Symbol}}
end

"""
    $(SIGNATURES)

EXPERIMENTAL: initialize target variable `xi` based on connected factors in the
factor graph `fgl`.  Possibly called from `addFactor!`, or `doCliqAutoInitUp!` (?).

Notes:
- Special carve out for multihypo cases, see issue 427.

Development Notes:
- Target factor is first (singletons) or second (dim 2 pairwise) variable vertex in `xi`.
- TODO use DFG properly with local operations and DB update at end.
- TODO get faster version of `isInitialized` for database version.
- TODO: Persist this back if we want to here.
- TODO: init from just partials
"""
function doautoinit!(
  dfg::AbstractDFG,
  xi::VariableCompute;
  solveKey::Symbol = :default,
  singles::Bool = true,
  N::Int = getSolverParams(dfg).N, #maximum([length(getPoints(getBelief(xi, solveKey))); getSolverParams(dfg).N]),
  logger = ConsoleLogger(),
)
  #
  didinit = false
  # create State if it does not exist yet
  prepareState!(xi, NPBPSolver(), solveKey; num_kernels=N)
  # don't initialize a variable more than once
  if !isInitialized(xi, solveKey)
    with_logger(logger) do
      @info "try doautoinit! of $(xi.label)"
    end
    # get factors attached to this variable xi
    vsym = xi.label
    neinodes = listNeighbors(dfg, vsym)
    # proceed if has more than one neighbor OR even if single factor
    if (singles || length(neinodes) > 1)
      # Which of the factors can be used for initialization
      useinitfct = Symbol[]
      # Consider factors connected to $vsym...
      for xifct in neinodes
        canuse, usefct, notusevars =
          factorCanInitFromOtherVars(dfg, xifct, vsym; solveKey = solveKey)
        if canuse
          union!(useinitfct, usefct)
        end
      end
      with_logger(logger) do
        @info "init with useinitfct $useinitfct"
      end
      # println("Consider all singleton (unary) factors to $vsym...")
      # calculate the predicted belief over $vsym
      if length(useinitfct) > 0
        with_logger(logger) do
          @info "do init of $vsym"
        end
        # FIXME ensure a product of only partial densities and returned pts are put to proper dimensions
        fcts = map(fx -> getFactor(dfg, fx), useinitfct)
        bel, ipc = propagateBelief(dfg, getVariable(dfg, vsym), fcts; solveKey, logger, N)
        # while the propagate step might allow large point counts, the graph should stay restricted to N
        bel_ =
          Npts(bel) == getSolverParams(dfg).N ? bel : resample(bel, getSolverParams(dfg).N)
        setBelief!(xi, bel_; solveKey) # TODO, update to stateLabel
        state = getState(xi, solveKey)
        state.initialized = true

        # Update the data in the event that it's not local
        # TODO perhaps use merge, but keeping to deepcopy as update variant used was set to copy.
        DFG.copytoState!(dfg, xi.label, solveKey, state)
        # deepcopy graphinit value, see IIF #612
        DFG.copytoState!(
          dfg,
          xi.label,
          :graphinit,
          getState(xi, solveKey),
        )
        didinit = true
      end
    end
  end
  return didinit
end

function doautoinit!(
  dfg::AbstractDFG,
  Xi::Vector{<:VariableCompute};
  solveKey::Symbol = :default,
  singles::Bool = true,
  N::Int = getSolverParams(dfg).N,
  logger = ConsoleLogger(),
)
  #
  #
  # Mighty inefficient function, since we only need very select fields nearby from a few neighboring nodes
  # do double depth search for variable nodes

  didinit = true

  # loop over all requested variables that must be initialized
  for xi in Xi
    didinit &=
      doautoinit!(dfg, xi; solveKey = solveKey, singles = singles, N = N, logger = logger)
  end
  return didinit
end

function doautoinit!(
  dfg::AbstractDFG,
  xsyms::Vector{Symbol};
  solveKey::Symbol = :default,
  singles::Bool = true,
  N::Int = getSolverParams(dfg).N,
  logger = ConsoleLogger(),
)
  #
  verts = getVariable.(dfg, xsyms)
  return doautoinit!(
    dfg,
    verts;
    solveKey = solveKey,
    singles = singles,
    N = N,
    logger = logger,
  )
end
function doautoinit!(
  dfg::AbstractDFG,
  xsym::Symbol;
  solveKey::Symbol = :default,
  singles::Bool = true,
  N::Int = getSolverParams(dfg).N,
  logger = ConsoleLogger(),
)
  #
  return doautoinit!(
    dfg,
    [getVariable(dfg, xsym);];
    solveKey = solveKey,
    singles = singles,
    N = N,
    logger = logger,
  )
end

"""
    $(TYPEDSIGNATURES)

Method to manually initialize a variable using a set of points.

Notes
- Disable automated graphinit on `addFactor!(fg, ...; graphinit=false)
  - any un-initialized variables will automatically be initialized by `solveTree!`

Example:

```julia
# some variable is added to fg
addVariable!(fg, :somepoint3, ContinuousEuclid{2})

# data is organized as (row,col) == (dimension, samples)
pts = randn(2,100)
initVariable!(fg, :somepoint3, pts)

# manifold management should be done automatically.
# note upgrades are coming to consolidate with Manifolds.jl, see RoME #244

## it is also possible to initVariable! by using existing factors, e.g.
initVariable!(fg, :x3, [:x2x3f1])
```

DevNotes
- TODO better document graphinit and treeinit.
"""
function initVariable!(
  variable::VariableCompute,
  ptsArr::ApproxManifoldProducts.HomotopyDensity,
  solveKey::Symbol = :default;
  # dontmargin::Bool = false,
  N::Int = length(getPoints(ptsArr)),
)
  #
  @debug "initVariable! $(getLabel(variable))"
  prepareState!(variable, NPBPSolver(), solveKey; num_kernels=N)
  setValKDE!(variable, ptsArr, true; solveKey)
  return nothing
end
function initVariable!(
  dfg::AbstractDFG,
  label::Symbol,
  belief::ApproxManifoldProducts.HomotopyDensity,
  solveKey::Symbol = :default;
  # dontmargin::Bool = false,
  N::Int = getSolverParams(dfg).N,
)
  #
  variable = getVariable(dfg, label)
  initVariable!(variable, belief, solveKey;
    # dontmargin = dontmargin,
    N = N
  )
  return nothing
end

function initVariable!(
  dfg::AbstractDFG,
  label::Symbol,
  samplable_belief::SamplableBelief,
  solveKey::Symbol = :default;
  N::Int = getSolverParams(dfg).N,
)
  #
  variable = getVariable(dfg, label)
  initVariable!(variable, samplable_belief, solveKey; N)
  return nothing
end

function initVariable!(
  variable::VariableCompute,
  samplable_belief::SamplableBelief,
  solveKey::Symbol = :default;
  N::Int = length(getVal(variable)),
)
  #
  M = getManifold(variable)
  if solveKey == :parametric
    prepareState!(variable, NLLSSolver(), solveKey)
    μ, iΣ = getMeasurementParametric(samplable_belief)
    vnd = getState(variable, solveKey)
    DFG.refMeans(vnd)[1] = getPoint(getStateKind(variable), μ)
    DFG.refCovariances(vnd)[1] .= inv(iΣ)
    vnd.initialized = true
  else
    points = [samplePoint(M, samplable_belief) for _ = 1:N]
    initVariable!(variable, points, solveKey)
  end
  return nothing
end

function initVariable!(
  dfg::AbstractDFG,
  label::Symbol,
  usefcts::AbstractVector{Symbol},
  solveKey::Symbol = :default;
  N::Int = getSolverParams(dfg).N,
  kwargs...,
)
  #
  pts = propagateBelief(dfg, label, usefcts; solveKey, N)[1]
  # pts = predictbelief(dfg, label, usefcts; solveKey = solveKey)[1]
  vert = getVariable(dfg, label)
  Xpre = manikde!(getManifold(getStateKind(vert)), pts)
  return initVariable!(vert, Xpre, solveKey; N, kwargs...)
  # setValKDE!(vert, Xpre, true, solveKey=solveKey)
  # return nothing
end

function initVariable!(
  vari::VariableCompute,
  pts::AbstractVector{P},
  solveKey::Symbol = :default;
  bw = nothing,
) where {P}
  #
  # specializations to support generic case of Tuple rather than ProductRepr or ArrayPartition inputs
  # TODO ArrayPartition inputs
  _prodrepr(pt) = pt
  # _prodrepr(pt::Tuple) = Manifolds.ProductRepr(pt...)
  _prodrepr(pt::Tuple) = ArrayPartition(pt...)
  _bw = bw === nothing ? zeros(getDimension(vari)) : bw
  M = getStateKind(vari)
  pp = HomotopyDensity_legacy(M, _prodrepr.(pts); bw=_bw, newbw=false)
  return initVariable!(vari, pp, solveKey)
end

function initVariable!(
  dfg::AbstractDFG,
  sym::Symbol,
  pts::AbstractVector{P},
  solveKey::Symbol = :default;
  kwargs...,
) where {P}
  #
  return initVariable!(getVariable(dfg, sym), pts, solveKey; kwargs...)
end

# legacy alias
const initVariableManual! = initVariable!

"""
    $SIGNATURES

Set solveKey values of `dest::AbstractDFG` according to `initKey::Symbol=:graphinit` values.

Notes
- Some flexibility for using two DFGs and different key values, see Examples and code for details.
- Can also be specific with `varList::Vector{Symbol}`.
- Returns `dest` graph.
- Uses the supersolve mechanism.

Examples
```julia
resetInitialValues!(fg)
resetInitialValues!(fg1,fg2)  # into 1 from 2
resetInitialValues!(fg1,fg1,:myotherinit)  # use different init value into solveKey :default
resetInitialValues!(fg1,fg1,:graphinit, :mysolver) # not into solveKey=:default but :mysolver
resetInitialValues!(fg1, varList=[:x1;:l3])  # Specific variables only

# Into `fgNew` object, leaving `fg` untouched
fgNew = deepcopy(fg)
resetInitialValues!(fgNew,fg)
```

Related

initVariable!, graphinit (keyword)
"""
function resetInitialValues!(
  dest::AbstractDFG,
  src::AbstractDFG = dest,
  initKey::Symbol = :graphinit,
  solveKey::Symbol = :default;
  varList::AbstractVector{Symbol} = ls(dest),
)
  #
  for vs in varList
    vnd = getState(getVariable(src, vs), initKey)
    # guess we definitely want to use copy to preserve the initKey memory
    # updateVariableSolverData!(dest, vs, vnd, solveKey, true; warn_if_absent = false)
    DFG.copytoState!(
      dest,
      vs,
      solveKey,
      vnd,
    )    
  end
  return dest
end
const resetInitValues! = resetInitialValues!

"""
    $SIGNATURES

Ensure that no variables set as `solvable=1` are floating free without any connected `solvable=1` factors.  If any found, then set those 'free' variable's `solvable=solvableFallback` (default `0`).

Related

[`initAll!`](@ref)
"""
function ensureSolvable!(
  dfg::AbstractDFG;
  solvableTarget::Int = 1,
  solvableFallback::Int = 0,
)
  # workaround in case isolated variables occur
  solvVars = ls(dfg; whereSolvable = >=(solvableTarget))
  varHasFact = (x -> length(listNeighbors(dfg, x; whereSolvable = >=(solvableTarget))) == 0).(solvVars)
  blankVars = solvVars[findall(varHasFact)]
  if 0 < length(blankVars)
    @warn(
      "solveTree! dissallows solvable variables without any connected solvable factors -- forcing solvable=0 on $(blankVars)"
    )
    (x -> setSolvable!(dfg, x, solvableFallback)).(blankVars)
  end
  return blankVars
end

"""
    $SIGNATURES

Prepare the factor graph for nonparametric belief propagation solving.

Specifically:
- Ensures every solvable variable has a `State` for `solveKey` with `N` particles (uninitialized).
- Ensures every solvable factor has a built CCW (`solvercache`).

No variable values are computed — use [`initAll!`](@ref) for that.

See also: [`initAll!`](@ref), [`prepareFactorCache!`](@ref)
"""
function prepare!(
  dfg::AbstractDFG,
  solver::NPBPSolver,
  statelabel::Symbol;
  num_kernels::Int = solver.defaultNumKernels,
  whereSolvable = >=(1),
)
  # Ensure all variables have State for statelabel
  count = prepareStates!(dfg, solver, statelabel; num_kernels, whereSolvable)

  # Ensure all factor caches (CCW) are built
  for fl in listFactors(dfg; whereSolvable)
    fct = getFactor(dfg, fl)
    if !isassigned(fct.solvercache)
      prepareFactorCache!(dfg, fct)
    end
  end

  return count
end

"""
    $SIGNATURES

Perform `graphinit` over all variables with `solvable=1` (default).


See also: [`ensureSolvable!`](@ref), [`prepare!`](@ref), (EXPERIMENTAL 'treeinit')
"""
function initAll!(
  dfg::AbstractDFG,
  solveKey::Symbol = :default;
  _parametricInit::Bool = solveKey === :parametric,
  solvable::Int = 1,
  N::Int = _parametricInit ? 1 : getSolverParams(dfg).N,
)
  #
  syms = intersect(DFG.getAddHistory(dfg), ls(dfg; whereSolvable = >=(solvable)))

  # Prepare all solver caches
  if _parametricInit
    prepare!(dfg, NLLSSolver(), solveKey)
  else
    prepare!(dfg, NPBPSolver(), solveKey; num_kernels = N, whereSolvable = >=(solvable))
  end

  # do the init
  repeatCount = 0
  repeatFlag = true
  while repeatFlag
    repeatFlag = false
    repeatCount += 1
    if 10 < repeatCount
      @info "not able to initialize all variables via the factor graph, abort autoinit."
      break
    end
    for sym in syms
      var = getVariable(dfg, sym)
      # is this SolverData initialized?
      if !isInitialized(var, solveKey)
        @info "$(var.label) is not initialized, and will do so now..."
        if _parametricInit
          autoinitParametric!(dfg, var; solveKey)
        else
          doautoinit!(dfg, [var;]; N, solveKey, singles = true)
        end
        !isInitialized(var, solveKey) ? (repeatFlag = true) : nothing
      end
    end
  end
  return nothing
end
