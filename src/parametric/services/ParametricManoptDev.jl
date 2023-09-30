using Manopt
using FiniteDiff
using SparseDiffTools
using BlockArrays
using SparseArrays

# using ForwardDiff
# using Zygote

##



function CalcFactorManopt(fg, fct::DFGFactor, varIntLabel)
  fac_func = getFactorType(fct)
  varOrder = getVariableOrder(fct)

  varOrderIdxs = getindex.(Ref(varIntLabel), varOrder)

  M = getManifold(getFactorType(fct))

  dims = manifold_dimension(M)

  meas, iΣ = getFactorMeasurementParametric(fct)

  sqrt_iΣ = convert(SMatrix{dims, dims}, sqrt(iΣ))
  cache = preambleCache(fg, getVariable.(fg, varOrder), getFactorType(fct))

  return CalcFactorManopt(
    fct.label,
    getFactorMechanics(fac_func),
    cache,
    varOrder,
    varOrderIdxs,
    meas,
    iΣ,
    sqrt_iΣ,
  )
end

function (cfm::CalcFactorManopt)(p)
  meas = cfm.meas
  idx = cfm.varOrderIdxs
  return cfm.sqrt_iΣ * cfm(meas, p[idx]...)
end

# cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
struct CostF_RLM!{T}
  points::Vector{T}
  costfuns::Vector{<:CalcFactorManopt}
end

function CostF_RLM!(costfuns::Vector{<:CalcFactorManopt}, frontals_p::Vector{T}, separators_p::Vector{T}) where T
  points::Vector{T} = vcat(frontals_p, separators_p)
  return CostF_RLM!(points, costfuns)
end

function (cfm::CostF_RLM!)(M::AbstractManifold, x::Vector, p::Vector{T}) where T
  cfm.points[1:length(p)] .= p
  return x .= mapreduce(f -> f(cfm.points), vcat, cfm.costfuns)
end

# jacobian of function for Riemannian Levenberg-Marquardt
struct JacF_RLM!{CF, T}
  costF!::CF
  X0::Vector{Float64}
  X::T
  q::T
  res::Vector{Float64}
end

# function JacF_RLM!(M, costF!; basis_domain::AbstractBasis = DefaultOrthonormalBasis())
function JacF_RLM!(M, costF!; basis_domain::AbstractBasis = DefaultOrthogonalBasis())
  
  p = costF!.points
  
  res = Vector(mapreduce(f -> f(p), vcat, costF!.costfuns))

  X0 = zeros(manifold_dimension(M))
  
  X = get_vector(M, p, X0, basis_domain)

  q = exp(M, p, X)

  # J = FiniteDiff.finite_difference_jacobian(
  #   Xc -> costF!(M, res, exp!(M, q, p, get_vector!(M, X, p, Xc, basis_domain))),
  #   X0,
  # )

  return JacF_RLM!(costF!, X0, X, q, res)

end

function (jacF!::JacF_RLM!)(
  M::AbstractManifold,
  J,
  p;
  # basis_domain::AbstractBasis = DefaultOrthonormalBasis(),
  basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
)
  
  X0 = jacF!.X0
  X = jacF!.X
  q = jacF!.q

  fill!(X0, 0)
  
    # ϵ = getPointIdentity(M)
    # function jaccost(res, Xc)
    #   exp!(M, q, ϵ, get_vector!(M, X, p, Xc, basis_domain))
    #   compose!(M, q, p, q)
    #   jacF!.costF!(M, res, q)
    # end


  # TODO maybe move to struct
  colorvec = matrix_colors(J)

  # TBD would the non-mutating staticarrays version be better?
  FiniteDiff.finite_difference_jacobian!(
    J,
    (res,Xc) -> jacF!.costF!(M, res, exp!(M, q, p, get_vector!(M, X, p, Xc, basis_domain))),
    X0;
    colorvec
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

function getSparsityPattern(fg)
  biadj = getBiadjacencyMatrix(fg)

  vdims = getDimension.(getVariable.(fg, biadj.varLabels))
  fdims = getDimension.(getFactor.(fg, biadj.facLabels))

  sm = map(eachindex(biadj.B)) do i
    vdim = vdims[i[2]]
    fdim = fdims[i[1]]
    if biadj.B[i] > 0
      trues(fdim,vdim)
    else
      falses(fdim,vdim)
    end
  end

  return SparseMatrixCSC(mortar(sm))

end

function solve_RLM(
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
  vartypes, vartypecount, vartypeslist = getVariableTypesCount(frontal_vars)
  
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
  fro_p = first.(getVal.(frontal_vars, solveKey = :parametric))
  sep_p::Vector{eltype(fro_p)} = first.(getVal.(fg, separators, solveKey = :parametric))

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostF_RLM!(calcfacs, fro_p, sep_p)
  # jacobian of function for Riemannian Levenberg-Marquardt
  jacF! = JacF_RLM!(MM, costF!)
  
  num_components = length(jacF!.res)

  p0 = deepcopy(fro_p)

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
    p0,
    num_components;
    evaluation=InplaceEvaluation(),
    initial_residual_values,
    initial_jacobian_f,
    kwargs...
  )

  return vartypeslist, lm_r
end

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
  jacF! = JacF_RLM!(MM, costF!)
  
  num_components = length(jacF!.res)

  p0 = deepcopy(fro_p)

  initial_residual_values = zeros(num_components)
  initial_jacobian_f = Float64.(getSparsityPattern(fg)) 

  #HEX solve
  # sparse J 0.025235 seconds (133.65 k allocations: 9.964 MiB
  # dense  J 0.022079 seconds (283.54 k allocations: 18.146 MiB)

  # 9.125818 seconds (86.35 M allocations: 6.412 GiB, 14.34% gc time)
  # 0.841720 seconds (7.96 M allocations: 751.825 MiB)

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
)
  @showprogress for vIdx in varorderIds
    autoinitParametricManopt!(fg, vIdx; reinit)
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

    vartypeslist, lm_r = solve_RLM(dfg, [initme], initfrom)

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
