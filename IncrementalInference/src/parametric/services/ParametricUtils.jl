# ================================================================================================
# Parametric Factors
# ================================================================================================

"""
    $SIGNATURES

Returns the parametric measurement for a factor as a tuple (measurement, inverse covariance) for parametric inference (assuming Gaussian).
Defaults to find the parametric measurement at field `Z`.
Notes
- Users should overload this method should their factor not default to `.Z<:ParametricType`.
- First design choice was to restrict this function to returning coordinates
  - See https://github.com/JuliaRobotics/RoME.jl/issues/465
  - Pay attention to which tangent space point is used for converting points on a manifold to coordinates,
    - Originally written just for Lie Groups to support legacy, but future needs may well alter the design.
- Original design driven by parametric solve and dead reckon tethering.

See also: [`accumulateFactorMeans`](@ref), [`solveFactorParametric`](@ref)
"""
function getMeasurementParametric end

function getMeasurementParametric(Z)
  return error(
    "$(typeof(Z)) is not supported, please use non-parametric or open an issue if it should be",
  )
end

function getMeasurementParametric(Z::Normal)
  meas = mean(Z)
  iσ = 1 / std(Z)^2
  return [meas], reshape([iσ], 1, 1)
end

function getMeasurementParametric(Z::MvNormal)
  meas = mean(Z)
  iΣ = invcov(Z)
  return meas, iΣ
end

# the point `p` on the manifold is the mean
function getMeasurementParametric(s::ManifoldPrior)
  meas = s.p
  iΣ = invcov(s.Z)
  return meas, iΣ
end

function getMeasurementParametric(s::AbstractObservation)
  if hasfield(typeof(s), :Z)
    Z = s.Z
  else
    error(
      "getMeasurementParametric(::$(typeof(s))) not defined, please add it, or use non-parametric, or open an issue for help.",
    )
  end

  return getMeasurementParametric(Z)
end

getMeasurementParametric(fct::FactorCompute) = getMeasurementParametric(getObservation(fct))
getMeasurementParametric(dfg::AbstractDFG, flb::Symbol) = getMeasurementParametric(getFactor(dfg, flb))

# maybe rename getMeasurementParametric to something like getNormalDistributionParams or getMeanCov

# default to point on manifold
function getFactorMeasurementParametric(fac::AbstractPriorObservation)
  M = getManifold(fac)
  ϵ = getPointIdentity(M)
  dims = manifold_dimension(M)
  Xc, iΣ = getMeasurementParametric(fac)
  # X = get_vector(M, ϵ, Xc, DefaultOrthogonalBasis())
  X = hat(LieAlgebra(M), Xc, typeof(ϵ))
  meas = convert(typeof(ϵ), exp(M, ϵ, X))
  iΣ = convert(SMatrix{dims, dims}, iΣ)
  meas, iΣ
end
# default to point on tangent vector
function getFactorMeasurementParametric(fac::AbstractRelativeObservation)
  M = getManifold(fac)
  ϵ = getPointIdentity(M)
  dims = manifold_dimension(M)
  Xc, iΣ = getMeasurementParametric(fac)
  # measX = convert(typeof(ϵ), get_vector(M, ϵ, Xc, DefaultOrthogonalBasis()))
  if M isa LieGroups.ValidationLieGroup
    measX = LieGroups.unwrap_validation(hat(LieAlgebra(M), Xc, typeof(ϵ)))
  else
    measX = convert(typeof(ϵ), hat(LieAlgebra(M), Xc, typeof(ϵ)))
  end
  iΣ = convert(SMatrix{dims, dims}, iΣ)
  measX, iΣ
end

getFactorMeasurementParametric(fct::FactorCompute) = getFactorMeasurementParametric(getObservation(fct))
getFactorMeasurementParametric(dfg::AbstractDFG, flb::Symbol) = getFactorMeasurementParametric(getFactor(dfg, flb))

function getFactorTypesCount(facs::Vector{<:FactorCompute})
  typedict = OrderedDict{DataType, Int}()
  alltypes = OrderedDict{DataType, Vector{Symbol}}()
  for f in facs
    facType = typeof(getObservation(f))
    cnt = get!(typedict, facType, 0)
    typedict[facType] = cnt + 1

    dt = get!(alltypes, facType, Symbol[])
    push!(dt, f.label)
  end
  #TODO tuple or vector?
  # vartypes = tuple(keys(typedict)...)
  factypes::Vector{DataType} = collect(keys(typedict))
  return factypes, typedict, alltypes
end

# ================================================================================================
# ================================================================================================
# New Parametric refactor WIP
# ================================================================================================
# ================================================================================================

# ================================================================================================
# GraphSolveStructures
# ================================================================================================

getVariableTypesCount(fg::AbstractDFG) = getVariableTypesCount(getVariables(fg))

function getVariableTypesCount(vars::Vector{<:VariableCompute})
  typedict = OrderedDict{DataType, Int}()
  alltypes = OrderedDict{DataType, Vector{Symbol}}()
  for v in vars
    varType = typeof(getStateKind(v))
    cnt = get!(typedict, varType, 0)
    typedict[varType] = cnt + 1

    dt = get!(alltypes, varType, Symbol[])
    push!(dt, v.label)
  end
  #TODO tuple or vector?
  # vartypes = tuple(keys(typedict)...)
  vartypes::Vector{DataType} = collect(keys(typedict))
  return vartypes, typedict, alltypes
end

buildGraphSolveManifold(fg::AbstractDFG) = buildGraphSolveManifold(getVariables(fg))

function buildGraphSolveManifold(vars::Vector{<:VariableCompute})
  vartypes, vartypecount, vartypeslist = getVariableTypesCount(vars)

  PMs = map(vartypes) do vartype
    N = vartypecount[vartype]
    G = getManifold(vartype)
    if G isa LieGroups.ValidationLieGroup
      #strip away ValidationLieGroup
      G = G.lie_group
    end
    return NPowerManifold(G, N)
    # return LieGroups.PowerLieGroup(G, N)
    # return Manifolds.PowerManifold(G, N)
    # PowerManifold(G, NestedReplacingPowerRepresentation(), N)
    # PowerManifold(G, NestedPowerRepresentation(), N) #TODO investigate as it does not converge
  end
  M = ProductManifold(PMs...)
  # M = ProductLieGroup(PMs...)
  return M, vartypes, vartypeslist
end

"""
    $SIGNATURES

Product manifold and label partition for variables taken in **partition order**: each group is grouped
by state type on its own, and the groups are concatenated.

Returns `(M, varlabelsAP, varIntLabel)`; the block order of `M` and of `varlabelsAP` agree, which is
the invariant that keeps a Jacobian's columns aligned with `Λ`'s coordinates.
"""
function buildPartitionedSolveManifold(varsets)
  PMs = []
  blocks = Vector{Symbol}[]
  for vars in varsets
    isempty(vars) && continue
    vartypes, vartypecount, vartypeslist = getVariableTypesCount(vars)
    for vartype in vartypes
      G = getManifold(vartype)
      if G isa LieGroups.ValidationLieGroup
        G = G.lie_group
      end
      push!(PMs, NPowerManifold(G, vartypecount[vartype]))
    end
    _, ap = getVarIntLabelMap(vartypeslist)
    append!(blocks, collect(ap.x))
  end
  M = ProductManifold(PMs...)
  blocks_tuple = (blocks...,)
  varlabelsAP = ArrayPartition{Symbol, typeof(blocks_tuple)}(blocks_tuple)
  varIntLabel = OrderedDict{Symbol, Int}(l => i for (i, l) in enumerate(varlabelsAP))
  return M, varlabelsAP, varIntLabel
end


function _get_dim_ranges(dims::NTuple{N,Any}) where {N}
  dims_acc = accumulate(+, vcat(1, SVector(dims)))
  return ntuple(i -> (dims_acc[i]:(dims_acc[i] + dims[i] - 1)), Val(N))
end

#NOTE this only works with a product of power manifolds
function getComponentsCovar(@nospecialize(PM::ProductManifold), Σ::AbstractMatrix)
  dims = manifold_dimension.(PM.manifolds)
  dim_ranges = _get_dim_ranges(dims)

  subsigmas = map(zip(dim_ranges, PM.manifolds)) do v
    r = v[1]
    M = v[2]
    return _getComponentsCovar(M, view(Σ, r, r))
  end

  return ArrayPartition(subsigmas...)
end

function _getComponentsCovar(@nospecialize(PM::PowerManifold), Σ::AbstractMatrix)
  M = PM.manifold
  dim = manifold_dimension(M)
  subsigmas = map(Manifolds.get_iterator(PM)) do i
    r = ((i - 1) * dim + 1):(i * dim)
    return Σ[r, r]
  end

  return subsigmas
end

function _getComponentsCovar(@nospecialize(PM::NPowerManifold), Σ::AbstractMatrix)
  M = PM.manifold
  dim = manifold_dimension(M)
  subsigmas = map(Manifolds.get_iterator(PM)) do i
    r = ((i - 1) * dim + 1):(i * dim)
    return Σ[r, r]
  end

  return subsigmas
end




# ================================================================================================
# Parametric utils
# ================================================================================================

# SANDBOX of usefull development functions to be cleaned up
"""
    $SIGNATURES
Update the parametric solver data value and covariance.
"""
function updateSolverDataParametric! end

function updateSolverDataParametric!(
  vnd::State,
  val::AbstractArray,
  cov::AbstractMatrix,
)
  # fill in the variable node data value
  DFG.refMeans(vnd)[1] = val
  #calculate and fill in covariance
  DFG.refCovariances(vnd)[1] .= cov
  return vnd
end

function updateSolverDataParametric!(
  v::VariableCompute,
  val::AbstractArray,
  cov::AbstractMatrix;
  solveKey::Symbol = :parametric,
)
  vnd = getState(v, solveKey)
  return updateSolverDataParametric!(vnd, val, cov)
end


"""
    $SIGNATURES
Initialize the parametric solver data from a different solution in `fromkey`.

DevNotes
- TODO, keyword `force` not wired up yet.
"""
function initParametricFrom!(
  fg::AbstractDFG,
  fromkey::Symbol = :default;
  parkey::Symbol = :parametric,
  onepoint = false,
  force::Bool = false,
)
  #
  if onepoint
    for v in getVariables(fg)
      fromvnd = getState(v, fromkey)
      dims = getDimension(v)
      DFG.refMeans(getState(v, parkey))[1] = DFG.refMeans(fromvnd)[1]
      DFG.refCovariances(getState(v, parkey))[1] = LinearAlgebra.I(dims)
    end
  else
    for var in getVariables(fg)
      dims = getDimension(var)
      μ, Σ = calcMeanCovar(var, fromkey)
      DFG.refMeans(getState(var, parkey))[1] = μ
      DFG.refCovariances(getState(var, parkey))[1] = Σ
    end
  end
end

"""
    $SIGNATURES
Add a parametric state label to all the variables in fg if it doesn't exist.
"""
function addParametricSolver!(fg; init = true, solveKey::Symbol = :parametric)
  if !(solveKey in fg.solverParams.algorithms)
    push!(fg.solverParams.algorithms, solveKey)
    foreach(
      v -> IIF.setDefaultNodeDataParametric!(v, getStateKind(v); solveKey, initialized = false),
      getVariables(fg),
    )
    if init
      autoinitParametric!(fg; solveKey)
    end
  else
    error("parametric solvekey $solveKey already exists")
  end
  return nothing
end

"""
    $SIGNATURES
Update the fg from solution in vardict. Usefull for plotting
"""
function updateParametricSolution!(sfg, vardict::AbstractDict; solveKey::Symbol = :parametric)
  for (v, val) in vardict
    vnd = getState(getVariable(sfg, v), solveKey)
    # Update the variable node data value and covariance
    updateSolverDataParametric!(vnd, val.val, val.cov)
  end
end

function updateParametricSolution!(fg, M, labels::AbstractArray{Symbol}, vals, Λ; solveKey::Symbol = :parametric)
  
  if isnothing(Λ)
    covars = nothing
  else
    Σ = try 
      cholesky(Symmetric(Matrix(Λ))) \ I
    catch ex
      if size(Λ, 1) < 1000 
        @warn "Precision matrix inversion failed, using pinv" ex
        pinv(Matrix(Λ))
      else
        @error "Precision matrix inversion failed and matrix is large, not updating covariance" ex
        nothing
      end
    end
    covars = isnothing(Σ) ? nothing : getComponentsCovar(M, Σ)
  end

  for (i, (v, val)) in enumerate(zip(labels, vals))
    vnd = getState(getVariable(fg, v), solveKey)
    covar = isnothing(covars) ? DFG.refCovariances(vnd)[1] : covars[i]
    # Update the variable node data value and covariance
    updateSolverDataParametric!(vnd, val, covar)
  end

end

function createMvNormal(val, cov)
  #TODO do something better for properly formed covariance, but for now just a hack...FIXME
  if all(diag(cov) .> 0.001) && isapprox(cov, transpose(cov); rtol = 1e-4)
    return MvNormal(val, Symmetric(cov))
  else
    @error("Covariance matrix error", cov)
    # return nothing # FIXME, blanking nothing during #459 consolidation
    return MvNormal(val, ones(length(val)))
  end
end

function createMvNormal(
  v::VariableCompute,
  key = :parametric;
  parametric::Bool = key === :parametric,
)
  if parametric
    vnd = getState(v, key)
    dims = getDimension(vnd)
    val = DFG.refMeans(vnd)[1]
    cov = DFG.refCovariances(vnd)[1:dims, 1:dims]
    return createMvNormal(val, cov)
  else
    @warn "Trying MvNormal Fit"
    return fit(MvNormal, DFG.refPoints(getState(v, key)))
  end
end

#TODO this is still experimental and a POC
"""
    $SIGNATURES

Build the Bayes tree for `fg` and return a vector of `(frontals, separators)` tuples
ordered root-to-leaves (BFS).
"""
function getInitOrderParametric(fg; ordering::Symbol = :qr)
  tree = buildTreeReset!(fg; ordering)

  # BFS root-to-leaves: parents are processed before children
  clique_order = Vector{NamedTuple{(:frontals, :separators), Tuple{Vector{Symbol}, Vector{Symbol}}}}()
  
  # Find root cliques
  queue = TreeClique[]
  for cliqId in getCliqueIds(tree)
    if isRoot(tree, cliqId)
      push!(queue, getClique(tree, cliqId))
    end
  end

  while !isempty(queue)
    cliq = popfirst!(queue)
    frontals = getCliqFrontalVarIds(cliq)
    separators = getCliqSeparatorVarIds(cliq)
    push!(clique_order, (; frontals, separators))
    # Enqueue children
    for child in getChildren(tree, cliq)
      push!(queue, child)
    end
  end

  return clique_order
end


#
