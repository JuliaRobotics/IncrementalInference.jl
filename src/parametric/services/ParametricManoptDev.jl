using Manopt
using FiniteDiff
using SparseDiffTools
using BlockArrays
using SparseArrays

# using ForwardDiff
# using Zygote

##
function getVarIntLabelMap(vartypeslist::OrderedDict{DataType, Vector{Symbol}})
  varlist_tuple = (values(vartypeslist)...,)
  varlabelsAP = ArrayPartition{Symbol, typeof(varlist_tuple)}(varlist_tuple)
  varIntLabel = OrderedDict(zip(varlabelsAP, collect(1:length(varlabelsAP))))
  return varIntLabel, varlabelsAP
end

function CalcFactorManopt(fg, fct::DFGFactor, varIntLabel, points::Union{Nothing,ArrayPartition}=nothing)
  fac_func = getFactorType(fct)
  varOrder = getVariableOrder(fct)

  varOrderIdxs = getindex.(Ref(varIntLabel), varOrder)

  M = getManifold(getFactorType(fct))

  dims = manifold_dimension(M)

  meas, iΣ = getFactorMeasurementParametric(fct)

  sqrt_iΣ = convert(SMatrix{dims, dims}, sqrt(iΣ))
  cache = preambleCache(fg, getVariable.(fg, varOrder), getFactorType(fct))

  if isnothing(points)
    points_view = nothing
  else
    points_view = @view points[varOrderIdxs]
  end

  return CalcFactorManopt(
    fct.label,
    getFactorMechanics(fac_func),
    cache,
    tuple(varOrder...),
    tuple(varOrderIdxs...),
    points_view,
    meas,
    iΣ,
    sqrt_iΣ,
  )
end

function calcFactorManoptVec(fg, factorLabels::Vector{Symbol}, varIntLabel::OrderedDict{Symbol, Int64}, points)
  factypes, typedict, alltypes = getFactorTypesCount(getFactor.(fg, factorLabels))
  
  # skip non-numeric prior (MetaPrior)
  #TODO test... remove MetaPrior{T} something like this
  metaPriorKeys = filter(k->contains(string(k), "MetaPrior"), collect(keys(alltypes)))
  delete!.(Ref(alltypes), metaPriorKeys)

  parts = map(values(alltypes)) do labels
    map(getFactor.(fg, labels)) do fct
      CalcFactorManopt(fg, fct, varIntLabel, points)
    end
  end
  parts_tuple = (parts...,)
  return ArrayPartition{CalcFactorManopt, typeof(parts_tuple)}(parts_tuple)
end

function (cfm::CalcFactorManopt{T})(p) where T
  meas = cfm.meas
  points = map(idx->p[idx], cfm.varOrderIdxs)
  return cfm.sqrt_iΣ * cfm(meas, points...) # 0.654783 seconds (6.75 M allocations: 531.688 MiB, 14.41% gc time)
end

function (cfm::CalcFactorManopt{T})() where T
  return cfm.sqrt_iΣ * cfm(cfm.meas, cfm.points...) # 0.654783 seconds (6.75 M allocations: 531.688 MiB, 14.41% gc time)
end

# cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
struct CostF_RLM!{T}
  points::Vector{T}
  costfuns::Vector{<:CalcFactorManopt}
  # varLabels::Vector{Symbol} @TODO add
end

function CostF_RLM!(costfuns::Vector{<:CalcFactorManopt}, frontals_p::Vector{T}, separators_p::Vector{T}) where T
  points::Vector{T} = vcat(frontals_p, separators_p)
  return CostF_RLM!(points, costfuns)
end

function (cfm::CostF_RLM!)(M::AbstractManifold, x::Vector, p::Vector) 
  cfm.points[1:length(p)] .= p
  # return x .= mapreduce(f -> f(cfm.points), vcat, cfm.costfuns)
  # return x .= reduce(vcat, map(f -> f(cfm.points), cfm.costfuns))
  st = 1
  for f in cfm.costfuns
    l = size(f.iΣ, 1)
    x[st:st + l - 1] = f(cfm.points)
    st += l
  end
end

struct CostF_RLM_AP!{CFT}
  # points::PT #TODO RENAME - don't update this in functor, seperator static points only!
  costfuns::ArrayPartition{CalcFactorManopt, CFT}
  varLabels::Vector{Symbol} # vector for performance above ArrayPartition{Symbol}?
  # varPoints::VPT
  # sepLabels::Vector{Symbol}
  # sepPoints::SPT
  # facLabels::Vector{Symbol}
  # add return_ranges to allow MultiThreaded
end

function calcFactorResVec!(x::Vector{T}, cfm_part::Vector{<:CalcFactorManopt}, p::AbstractArray{T}, st::Int) where T
  l = getDimension(cfm_part[1]) # all should be the same 
  for cfm in cfm_part
    x[st:st + l - 1] = cfm(p) #NOTE looks like do not broadcast here
    st += l
  end
  return st
end

function calcFactorResVec_threaded!(x::Vector{T}, cfm_part::Vector{<:CalcFactorManopt}, p::AbstractArray{T}, st::Int) where T
  l = getDimension(cfm_part[1]) # all should be the same 
  Threads.@threads for i in eachindex(cfm_part)
    r = range(st+l*(i-1), length=l)
    cfm = cfm_part[i]
    x[r] = cfm(p) #NOTE looks like do not broadcast here
  end
  return st + l*length(cfm_part)
end

# function (cfm::CostF_RLM!)(M::AbstractManifold, x::Vector, p::Vector{T}) where T
function (costf::CostF_RLM_AP!{CFT})(M::AbstractManifold, x::Vector{T}, p::AbstractVector{T}) where {CFT,T}
  # costf.points[1:length(p)] .= p #FIXME remove just testing sparse kwarg
  st = 1
  for cfm_part in costf.costfuns.x
      st = calcFactorResVec!(x, cfm_part, p, st)
  end
  return x
end

## --------------------------------------------------------------------------------------------------------------
## jacobian of function for Riemannian Levenberg-Marquardt
## --------------------------------------------------------------------------------------------------------------
struct JacF_RLM!{CF, T, JC}
  costF!::CF
  X0::Vector{Float64}
  X::T
  q::T
  res::Vector{Float64}
  Jcache::JC
end

# function JacF_RLM!(M, costF!; basis_domain::AbstractBasis = DefaultOrthonormalBasis())
function JacF_RLM!(M, costF!, p, fg=nothing; basis_domain::AbstractBasis = DefaultOrthogonalBasis(), sparse=!isnothing(fg))
  
  #Why does this error?
  # res = Vector(mapreduce(f -> f(p), vcat, costF!.costfuns)) 
  res = reduce(vcat, map(f -> f(p), Vector(costF!.costfuns)))

  X0 = zeros(manifold_dimension(M))
  
  X = get_vector(M, p, X0, basis_domain)

  q = exp(M, p, X)

  if sparse
    factLabels = collect(getproperty.(costF!.costfuns, :faclbl))
    sparsity = eltype(res).(getSparsityPattern(fg, costF!.varLabels, factLabels))
    colorvec = matrix_colors(sparsity)
  else 
    sparsity = nothing
    colorvec = 1:length(X0)
  end

  cache = FiniteDiff.JacobianCache(X0, res; colorvec, sparsity)

  return JacF_RLM!(costF!, X0, X, q, res, cache)

end

# TODO addd M to JacF_RLM! and test this ipo closure
# function (jacF!::JacF_RLM!)(res, Xc)
#   X = jacF!.X
#   q = jacF!.q
#   get_vector!(M, X, p, Xc, basis_domain)
#   exp!(M, q, p, X)
#   return jacF!.costF!(M, res, q)
# end

function (jacF!::JacF_RLM!)(
  M::AbstractManifold,
  J,
  p::T;
  # basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
  basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
) where T
  
  X0 = jacF!.X0
  X = jacF!.X
  q = jacF!.q
  cache = jacF!.Jcache
  
  fill!(X0, 0)
  
    # ϵ = getPointIdentity(M)
    # function jaccost(res, Xc)
    #   exp!(M, q, ϵ, get_vector!(M, X, p, Xc, basis_domain))
    #   compose!(M, q, p, q)
    #   jacF!.costF!(M, res, q)
    # end

  # TODO make sure closure performs (let, ::, or (jacF!::JacF_RLM!)(res, Xc))
  function costf!(res, Xc)
    get_vector!(M, X, p, Xc, basis_domain)
    exp!(M, q, p, X)
    jacF!.costF!(M, res, q)
  end

  # colorvec = matrix_colors(J) #FIXME, this should work from cache, but id doesn't 

  FiniteDiff.finite_difference_jacobian!(
    J,
    costf!,
    X0,
    cache;
    # colorvec #FIXME also should work from cache
  )

  # ManifoldDiff._jacobian!(
  #   J, 
  #   (Xc)->jacF!.costF!(M, jacF!.res, exp!(M, q, p, get_vector!(M, X, p, Xc, basis_domain))),
  #   X0,
  #   ManifoldDiff.default_differential_backend()
  # )

  return J
end

struct FactorGradient{A <: AbstractMatrix}
  manifold::AbstractManifold
  JacF!::JacF_RLM!
  J::A
end

# TODO this function is not the sparsity pattern yet, it just fills in all entries from the biadjacency matrix
# TODO allow getting sparcity pattern for a subfg
# OLD 0.424040 seconds (940.11 k allocations: 45.512 MiB)
# NEW 0.001552 seconds (2.04 k allocations: 1.816 MiB)
function getSparsityPattern(fg, varLabels, factLabels)
  biadj = getBiadjacencyMatrix(fg; varLabels, factLabels)

  vdims = getDimension.(getVariable.(fg, biadj.varLabels))
  fdims = getDimension.(getFactor.(fg, biadj.facLabels))
 
  c_end = cumsum(vdims)
  r_end = cumsum(fdims)

  C_range = range.(c_end - vdims .+1, c_end)
  R_range = range.(r_end - fdims .+1, r_end)

  ROWS, COLS, _ = findnz(biadj.B)

  iter = reduce(vcat, map(zip(ROWS, COLS)) do (R,C)
    vec(CartesianIndices((R_range[R], C_range[C])))
  end)

  vec(CartesianIndices((R_range[2], C_range[1])))

  return sparse(getindex.(iter,1), getindex.(iter,2), ones(Bool, length(iter)))
end


function solve_RLM(
  fg,
  frontals::Vector{Symbol} = ls(fg);
  sparse=false,
  # separators::Vector{Symbol} = setdiff(ls(fg), frontals);
  kwargs...
)
  #FIXME
  separators = Symbol[]
  # get the subgraph formed by all frontals, separators and fully connected factors
  varlabels = union(frontals, separators)
  faclabels = sortDFG(setdiff(getNeighborhood(fg, varlabels, 1), varlabels))

  filter!(faclabels) do fl
    return issubset(getVariableOrder(fg, fl), varlabels)
  end

  # so the subgraph consists of varlabels(frontals + separators) and faclabels

  # varIntLabel = OrderedDict(zip(varlabels, collect(1:length(varlabels))))

  # varIntLabel_frontals = filter(p->first(p) in frontals, varIntLabel)
  # varIntLabel_separators = filter(p->first(p) in separators, varIntLabel)
 
  # get the manifold and variable types
  frontal_vars = getVariable.(fg, frontals)
   
  M, varTypes, vartypeslist = buildGraphSolveManifold(frontal_vars)

  varIntLabel, varlabelsAP = getVarIntLabelMap(vartypeslist)

  #Can use varIntLabel (because its an OrderedDict), but varLabelsAP makes the ArrayPartition.
  # TODO just make sure this map makes a copy
  p0 = map(varlabelsAP) do label
    getVal(fg, label, solveKey = :parametric)[1]
  end

  if isempty(separators)

    all_points = p0
    
  else
    error("#FIXME this does not work yet with separators")
    _, _, all_vartypeslist = getVariableTypesCount(getVariable.(fg,varlabels))
    all_varIntLabel, varlabelsAP = getVarIntLabelMap(all_vartypeslist)

    all_points = map(varlabelsAP) do label
      getVal(fg, label, solveKey = :parametric)[1]
    end

  end

  #FIXME split frontals and separators
  # create an ArrayPartition{CalcFactorManopt} for faclabels
  calcfacs = calcFactorManoptVec(fg, faclabels, varIntLabel, all_points)

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostF_RLM_AP!(calcfacs, collect(varlabelsAP))

  @debug "building jacF!"
  # jacobian of function for Riemannian Levenberg-Marquardt
  jacF! = JacF_RLM!(M, costF!, p0, fg; sparse)

  num_components = length(jacF!.res)
  initial_residual_values = zeros(num_components)

  # initial_jacobian_f not type stable, but function barrier so should be ok.
  initial_jacobian_f = sparse ? 
    jacF!.Jcache.sparsity : 
    zeros(num_components, manifold_dimension(M))

  lm_r = Manopt.LevenbergMarquardt!(
    M,
    costF!,
    jacF!,
    p0,
    num_components;
    evaluation=InplaceEvaluation(),
    jacobian_tangent_basis = DefaultOrthogonalBasis(),
    initial_residual_values,
    initial_jacobian_f,
    kwargs...
  )

  return varIntLabel, lm_r
end

  # nlso = NonlinearLeastSquaresObjective(
  #   costF!,
  #   jacF!,
  #   num_components;
  #   evaluation = InplaceEvaluation(),
  #   jacobian_tangent_basis = DefaultOrthogonalBasis(),
  # )

  # @debug "starting solver"
  # lm_r = LevenbergMarquardt!(
  #   M, nlso, p0; 
  #   evaluation = InplaceEvaluation(),
  #   jacobian_tangent_basis = DefaultOrthogonalBasis(),
  #   initial_residual_values,
  #   initial_jacobian_f,
  #   kwargs...
  # )


#TODO implement frontals and seperators in new for init
function solve_RLM_conditional(
  fg,
  frontals::Vector{Symbol} = ls(fg),
  separators::Vector{Symbol} = setdiff(ls(fg), frontals);
  kwargs...
)

  # get the subgraph formed by all frontals, separators and fully connected factors
  varlabels = union(frontals, separators)
  faclabels = sortDFG(setdiff(getNeighborhood(fg, varlabels, 1), varlabels))

  filter!(faclabels) do fl
    return issubset(getVariableOrder(fg, fl), varlabels)
  end

  facs = getFactor.(fg, faclabels)

  # so the subgraph consists of varlabels(frontals + separators) and faclabels

  varIntLabel = OrderedDict(zip(varlabels, collect(1:length(varlabels))))

  # varIntLabel_frontals = filter(p->first(p) in frontals, varIntLabel)
  # varIntLabel_separators = filter(p->first(p) in separators, varIntLabel)

  calcfacs = map(f->CalcFactorManopt(fg, f, varIntLabel), facs)

  # get the manifold and variable types
  frontal_vars = getVariable.(fg, frontals)
   
  M, varTypes, vartypeslist = buildGraphSolveManifold(frontal_vars)
  #
  #FIXME 
  @assert length(M.manifolds) == 1 "#FIXME, this only works with 1 manifold type component"
  MM = M.manifolds[1]

  # inital values and separators from fg
  fro_p = first.(getVal.(frontal_vars, solveKey = :parametric))
  sep_p::Vector{eltype(fro_p)} = first.(getVal.(fg, separators, solveKey = :parametric))

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostF_RLM!(calcfacs, fro_p, sep_p)
  # jacobian of function for Riemannian Levenberg-Marquardt
  jacF! = JacF_RLM!(MM, costF!, vcat(fro_p, sep_p))
  
  num_components = length(jacF!.res)

  initial_residual_values = zeros(num_components)
  initial_jacobian_f = zeros(num_components, manifold_dimension(MM))
# 
  #HEX solve
  # sparse J 0.025235 seconds (133.65 k allocations: 9.964 MiB
  # dense  J 0.022079 seconds (283.54 k allocations: 18.146 MiB)

  lm_r = LevenbergMarquardt(
    MM,
    costF!,
    jacF!,
    fro_p,
    num_components;
    evaluation=InplaceEvaluation(),
    initial_residual_values,
    initial_jacobian_f,
    kwargs...
  )

  return vartypeslist, lm_r
end

  #HEX solve
  # sparse J 0.025235 seconds (133.65 k allocations: 9.964 MiB
  # new1     0.013486 seconds (36.16 k allocations: 2.593 MiB)
  # new2    0.010764 seconds (34.61 k allocations: 3.111 MiB)
  # dense  J 0.022079 seconds (283.54 k allocations: 18.146 MiB)
  
function solve_RLM_sparse(fg; kwargs...)

  # get the subgraph formed by all frontals, separators and fully connected factors
  varlabels = ls(fg)
  faclabels = lsf(fg)

  facs = getFactor.(fg, faclabels)

  # so the subgraph consists of varlabels(frontals + separators) and faclabels

  varIntLabel = OrderedDict(zip(varlabels, collect(1:length(varlabels))))

  calcfacs = CalcFactorManopt.(fg, facs, Ref(varIntLabel))

  # get the manifold and variable types
  vars = getVariable.(fg, varlabels)
  vartypes, vartypecount, vartypeslist = getVariableTypesCount(vars)
  
  PMs = map(vartypes) do vartype
    N = vartypecount[vartype]
    G = getManifold(vartype)
    return IIF.NPowerManifold(G, N)
  end
  M = ProductManifold(PMs...)
  
  #
  #FIXME 
  @assert length(M.manifolds) == 1 "#FIXME, this only works with 1 manifold type component"
  MM = M.manifolds[1]

  # inital values and separators from fg
  fro_p = first.(getVal.(vars, solveKey = :parametric))
  sep_p = eltype(fro_p)[]

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostF_RLM!(calcfacs, fro_p, sep_p)
  # jacobian of function for Riemannian Levenberg-Marquardt
  jacF! = JacF_RLM!(MM, costF!, fro_p, fg; sparse=true)
  
  num_components = length(jacF!.res)

  p0 = deepcopy(fro_p)

  initial_residual_values = zeros(num_components)
  initial_jacobian_f = jacF!.Jcache.sparsity 

  lm_r = LevenbergMarquardt(
    MM,
    costF!,
    jacF!,
    p0,
    num_components;
    evaluation=InplaceEvaluation(),
    initial_residual_values,
    initial_jacobian_f,
    kwargs...
  )

  return vartypeslist, lm_r
end

function autoinitParametricManopt!(
  fg,
  varorderIds = getInitOrderParametric(fg);
  reinit = false,
  kwargs...
)
  @showprogress for vIdx in varorderIds
    autoinitParametricManopt!(fg, vIdx; reinit, kwargs...)
  end
  return nothing
end

function autoinitParametricManopt!(dfg::AbstractDFG, initme::Symbol; kwargs...)
  return autoinitParametricManopt!(dfg, getVariable(dfg, initme); kwargs...)
end

function autoinitParametricManopt!(
  dfg::AbstractDFG,
  xi::DFGVariable;
  solveKey = :parametric,
  reinit::Bool = false,
  kwargs...,
)
  #

  initme = getLabel(xi)
  vnd = getSolverData(xi, solveKey)
  # don't initialize a variable more than once
  if reinit || !isInitialized(xi, solveKey)

    # frontals - initme
    # separators - inifrom

    initfrom = ls2(dfg, initme)
    filter!(initfrom) do vl
      return isInitialized(dfg, vl, solveKey)
    end

    vartypeslist, lm_r = solve_RLM_conditional(dfg, [initme], initfrom; kwargs...)

    val = lm_r[1]
    vnd::VariableNodeData = getSolverData(xi, solveKey)
    vnd.val[1] = val
  
    
    # val = lm_r[1]
    # cov =  ...
    # updateSolverDataParametric!(vnd, val, cov)

    vnd.initialized = true
    #fill in ppe as mean
    Xc::Vector{Float64} = collect(getCoordinates(getVariableType(xi), val))
    ppe = MeanMaxPPE(:parametric, Xc, Xc, Xc)
    getPPEDict(xi)[:parametric] = ppe

    result = vartypeslist, lm_r

  else
    result = nothing
  end

  return result#isInitialized(xi, solveKey)
end


## Check when time and delete if it can't be improved, curretnly ArrayPartition works best
#=
using FunctionWrappers: FunctionWrapper

# call with 
calcfacs = calcFactorManoptWrapper(fg, faclabels, varIntLabel, all_points)
costF! = CostF_RLM_WRAP2!(all_points, calcfacs, map(cfm->size(cfm.obj.x.iΣ,1), calcfacs))

function calcFactorManoptWrapper(fg, factorLabels::Vector{Symbol}, varIntLabel::OrderedDict{Symbol, Int64}, points::ArrayPartition)
  factypes, typedict, alltypes = getFactorTypesCount(getFactor.(fg, factorLabels))
  
  # skip non-numeric prior (MetaPrior)
  #TODO test... remove MetaPrior{T} something like this
  metaPriorKeys = filter(k->contains(string(k), "MetaPrior"), collect(keys(alltypes)))
  delete!.(Ref(alltypes), metaPriorKeys)

  calcfacs = map(factorLabels) do labels
    fct = getFactor(fg, labels)
    # moet ek 'n view in p0 in maak wat jy net p0 update en CFM view automaties, toets dit...
    cfm = IIF.CalcFactorManopt(fg, fct, varIntLabel, points)
    # return FunctionWrapper{Vector{Float64}, Tuple{typeof(points)}}(cfm)
    return FunctionWrapper{Vector{Float64}, Tuple{}}(cfm)
  end
  return calcfacs
end

struct CostF_RLM_WRAP2!{PT, CFW}
  points::PT
  costfuns::Vector{CFW}
  retdims::Vector{Int}
end

# function (cfm::CostF_RLM!)(M::AbstractManifold, x::Vector, p::Vector{T}) where T
function (cost::CostF_RLM_WRAP2!)(M::AbstractManifold, x::Vector{T}, p::AbstractVector{T}) where T
  # x .= reduce(vcat, map(f -> f(p), cost.costfuns))
  # x .= reduce(vcat, map(f -> f(), cost.costfuns))
  st = 1
  for (d, f) in zip(cost.retdims, cost.costfuns)
    x[st:st + d - 1] .= f(p)
    # x[st:st + d - 1] .= f()
    # fx = f.obj.x
    # x[st:st + d - 1] = fx.sqrt_iΣ * fx(fx.meas, fx.points...)
    st += d
  end
  return x
end
=#

