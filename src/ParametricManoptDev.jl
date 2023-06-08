using Manopt
using FiniteDiff
# using ForwardDiff
# using Zygote

##

struct CalcFactorManopt{
  D,
  L,
  FT <: AbstractFactor,
  M <: AbstractManifold,
  MEAS <: AbstractArray,
}
  faclbl::Symbol
  calcfactor!::CalcFactor{FT, Nothing, Nothing, Tuple{}, M}
  varOrder::Vector{Symbol}
  varOrderIdxs::Vector{Int}
  meas::MEAS
  iΣ::SMatrix{D, D, Float64, L}
end

function CalcFactorManopt(fg, fct::DFGFactor, varIntLabel)
  fac_func = getFactorType(fct)
  varOrder = getVariableOrder(fct)

  varOrderIdxs = getindex.(Ref(varIntLabel), varOrder)

  M = getManifold(getFactorType(fct))

  dims = manifold_dimension(M)
  ϵ = getPointIdentity(M)

  _meas, _iΣ = getMeasurementParametric(fct)
  if fac_func isa ManifoldPrior
    meas = _meas # already a point on M
  elseif fac_func isa AbstractPrior
    X = get_vector(M, ϵ, _meas, DefaultOrthogonalBasis())
    meas = exp(M, ϵ, X) # convert to point on M
  else
    # its a relative factor so should be a tangent vector 
    meas = convert(typeof(ϵ), get_vector(M, ϵ, _meas, DefaultOrthogonalBasis()))
  end

  # make sure its an SMatrix
  iΣ = convert(SMatrix{dims, dims}, _iΣ)

  # cache = preambleCache(fg, getVariable.(fg, varOrder), getFactorType(fct))

  calcf = CalcFactor(
    getFactorMechanics(fac_func),
    0,
    nothing,
    true,
    nothing,#cache,
    (), #DFGVariable[],
    0,
    getManifold(fac_func),
  )
  return CalcFactorManopt(fct.label, calcf, varOrder, varOrderIdxs, meas, iΣ)
end

function (cfm::CalcFactorManopt)(p)
  meas = cfm.meas
  idx = cfm.varOrderIdxs
  return cfm.calcfactor!(meas, p[idx]...)
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

function (cfm::CostF_RLM!)(M::AbstractManifold, x, p::Vector{T}) where T
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

function JacF_RLM!(M, costF!; basis_domain::AbstractBasis = DefaultOrthogonalBasis())
  
  p = costF!.points
  
  res = mapreduce(f -> f(p), vcat, costF!.costfuns)

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
  basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
)
  
  X0 = jacF!.X0
  X = jacF!.X
  q = jacF!.q

  fill!(X0, 0)

  # J .= FiniteDiff.finite_difference_jacobian(
  #   Xc -> jacF!.costF!(M, jacF!.res, exp!(M, q, p, get_vector!(M, X, p, Xc, basis_domain))),
  #   X0,
  # )
  FiniteDiff.finite_difference_jacobian!(
    J,
    (res,Xc) -> jacF!.costF!(M, res, exp!(M, q, p, get_vector!(M, X, p, Xc, basis_domain))),
    X0,
  )
  return J
end


function solve_RLM(
  fg,
  frontals::Vector{Symbol} = ls(fg),
  separators::Vector{Symbol} = setdiff(ls(fg), frontals);
)
  @error "#FIXME, use covariances" maxlog=1

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

  calcfacs = CalcFactorManopt.(fg, facs, Ref(varIntLabel))

  
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


  if false
    # non-in-place version updated bolow to in-place
    fullp::Vector{eltype(fro_p)} = [fro_p; sep_p] 
    # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
    function costF_RLM(M::AbstractManifold, p::Vector{T}) where T
      fullp[1:length(p)] .= p
      return Vector(mapreduce(f -> f(fullp), vcat, calcfacs))
    end

    # jacobian of function for Riemannian Levenberg-Marquardt
    function jacF_RLM(
      M::AbstractManifold,
      p;
      basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
    )
      X0 = zeros(manifold_dimension(M))
      J = FiniteDiff.finite_difference_jacobian(
        x -> costF_RLM(M, exp(M, p, get_vector(M, p, x, basis_domain))),
        X0,
      )
      # J = ForwardDiff.jacobian(
      #     x -> costF_RLM(M, exp(M, p, get_vector(M, p, x, basis_domain))),
      #     X0,
      # )
      return J
    end

    # 0.296639 seconds (2.46 M allocations: 164.722 MiB, 12.83% gc time)
    p0 = deepcopy(fro_p)
    lm_r = LevenbergMarquardt(MM, costF_RLM, jacF_RLM, p0)

    # 81.185117 seconds (647.20 M allocations: 41.680 GiB, 8.61% gc time)
  else
    # 74.420872 seconds (567.70 M allocations: 34.698 GiB, 8.30% gc time, 0.66% compilation time)
    #cost and jacobian functions
    # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
    costF! = CostF_RLM!(calcfacs, fro_p, sep_p)
    # jacobian of function for Riemannian Levenberg-Marquardt
    jacF! = JacF_RLM!(MM, costF!)
    
    num_components = length(jacF!.res)

    p0 = deepcopy(fro_p)
    lm_r = LevenbergMarquardt(MM, costF!, jacF!, p0, num_components; evaluation=InplaceEvaluation())
  end

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
