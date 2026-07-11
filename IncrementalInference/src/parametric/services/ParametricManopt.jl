using Manopt
using Manopt.Printf
using FiniteDiff
using SparseDiffTools
using SparseArrays

using ForwardDiff
# using Zygote

##
function getVarIntLabelMap(
  vartypeslist::OrderedDict{DataType, Vector{Symbol}}
)
  varlist_tuple = (values(vartypeslist)...,)
  varlabelsAP = ArrayPartition{Symbol, typeof(varlist_tuple)}(varlist_tuple)
  varIntLabel = OrderedDict(zip(varlabelsAP, collect(1:length(varlabelsAP))))
  return varIntLabel, varlabelsAP
end

function CalcFactorResidual(
  fg, 
  fct::FactorCompute, 
  varIntLabel
)
  fac_func = getObservation(fct)
  varOrder = collect(getVariableOrder(fct))

  varOrderIdxs = getindex.(Ref(varIntLabel), varOrder)

  M = getManifold(getObservation(fct))

  dims = manifold_dimension(M)

  meas, iΣ = getFactorMeasurementParametric(fct)

  sqrt_iΣ = convert(SMatrix{dims, dims}, sqrt(iΣ))
  cache = preambleCache(fg, getVariable.(fg, varOrder), getObservation(fct))

  return CalcFactorResidual(
    fct.label,
    fac_func,
    tuple(varOrder...),
    tuple(varOrderIdxs...),
    meas,
    sqrt_iΣ,
    cache,
  )
end


"""
  CalcFactorResidualAP
Create an `ArrayPartition` of `CalcFactorResidual`s.
"""
function CalcFactorResidualAP(
  fg::GraphsDFG, 
  factorLabels::Vector{Symbol}, 
  varIntLabel::OrderedDict{Symbol, Int64}
)
  factypes, typedict, alltypes = getFactorTypesCount(getFactor.(fg, factorLabels))
  
  # skip non-numeric prior (MetaPrior)
  #TODO test... remove MetaPrior{T} something like this
  metaPriorKeys = filter(k->contains(string(k), "MetaPrior"), collect(keys(alltypes)))
  delete!.(Ref(alltypes), metaPriorKeys)

  parts = map(values(alltypes)) do labels
    map(getFactor.(fg, labels)) do fct
      CalcFactorResidual(fg, fct, varIntLabel)
    end
  end
  parts_tuple = (parts...,)
  return ArrayPartition{CalcFactorResidual, typeof(parts_tuple)}(parts_tuple)
end

function (cfm::CalcFactorResidual)(p)
  meas = cfm.meas
  points = map(idx->p[idx], cfm.varOrderIdxs)
  return cfm.sqrt_iΣ * cfm(meas, points...)
end

# cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
struct CostFres_cond!{PT, CFT}
  points::PT
  costfuns::ArrayPartition{CalcFactorResidual, CFT}
  varLabels::Vector{Symbol}
end

function (costf::CostFres_cond!)(M::AbstractManifold, x::Vector, p::AbstractVector) 
  
  costf.points[1:length(p)] .= p

  st = 1
  for cfm_part in costf.costfuns.x
    st = calcFactorResVec!(x, cfm_part, costf.points, st)
  end
  return x

end

struct CostFres!{CFT}
  # points::PT #TODO RENAME - don't update this in functor, seperator static points only!
  costfuns::ArrayPartition{CalcFactorResidual, CFT}
  varLabels::Vector{Symbol} # vector for performance above ArrayPartition{Symbol}?
  # varPoints::VPT
  # sepLabels::Vector{Symbol}
  # sepPoints::SPT
  # facLabels::Vector{Symbol}
  # add return_ranges to allow MultiThreaded
end

function calcFactorResVec!(
  x::Vector{T},
  cfm_part::Vector{<:CalcFactorResidual{FT, N, D}},
  p::AbstractArray,
  st::Int
) where {T, FT, N, D}
  for cfm in cfm_part
    x[st:st + D - 1] = cfm(p) #NOTE looks like do not broadcast here
    st += D
  end
  return st
end

function calcFactorResVec_threaded!(x::Vector{T}, cfm_part::Vector{<:CalcFactorResidual}, p::AbstractArray, st::Int) where T
  l = getDimension(cfm_part[1]) # all should be the same
  N = length(cfm_part)
  chunkies = Iterators.partition(1:N, N ÷ Threads.nthreads())
  Threads.@threads for chunki in collect(chunkies)
    for i in chunki
      r = range(st + l*(i - 1); length = l)
      cfm = cfm_part[i]
      x[r] = cfm(p) #NOTE looks like do not broadcast here
    end
  end
  return st + l*N
end

function (costf::CostFres!{CFT})(M::AbstractManifold, x::Vector{T}, p::AbstractVector) where {CFT,T}
  st = 1
  for cfm_part in costf.costfuns.x
      # if length(cfm_part) > Threads.nthreads() * 10
        # st = calcFactorResVec_threaded!(x, cfm_part, p, st)
      # else
        st = calcFactorResVec!(x, cfm_part, p, st)
      # end
  end
  return x
end

## --------------------------------------------------------------------------------------------------------------
## jacobian of function for Riemannian Levenberg-Marquardt
## --------------------------------------------------------------------------------------------------------------
struct JacF_RLM!{CF, TX, TQ, JC}
  costF!::CF
  X0::Vector{Float64}
  X::TX
  q::TQ
  res::Vector{Float64}
  Jcache::JC
end

# function JacF_RLM!(M, costF!; basis_domain::AbstractBasis = DefaultOrthonormalBasis())
function JacF_RLM!(M, costF!, p, fg=nothing;
  all_points=p,
  basis_domain::AbstractBasis = LieGroups.DefaultLieAlgebraOrthogonalBasis(),
  is_sparse=!isnothing(fg)
)

  res = reduce(vcat, map(f -> f(all_points), Vector(costF!.costfuns)))

  X0 = zeros(manifold_dimension(M))
  
  X = get_vector(M, p, X0, basis_domain)
  # X = vee(LieAlgebra(M), X0)

  q = exp(M, p, X)

  if is_sparse
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
  # basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
  basis_domain::AbstractBasis = LieGroups.DefaultLieAlgebraOrthogonalBasis(),
) where T
  
  X0 = jacF!.X0
  X = jacF!.X
  q = jacF!.q
  cache = jacF!.Jcache
  
  fill!(X0, 0)

  # TODO make sure closure performs (let, ::, or (jacF!::JacF_RLM!)(res, Xc))
  function costf!(res, Xc)
    get_vector!(M, X, p, Xc, basis_domain)
    exp!(M, q, p, X)
    jacF!.costF!(M, res, q)
  end

  FiniteDiff.finite_difference_jacobian!(
    J,
    costf!,
    X0,
    cache;
  )

  return J
end

  # ϵ = getPointIdentity(M)
  # function jaccost(res, Xc)
  #   exp!(M, q, ϵ, get_vector!(M, X, p, Xc, basis_domain))
  #   compose!(M, q, p, q)
  #   jacF!.costF!(M, res, q)
  # end

  # ManifoldDiff._jacobian!(
  #   J, 
  #   (Xc)->jacF!.costF!(M, jacF!.res, exp!(M, q, p, get_vector!(M, X, p, Xc, basis_domain))),
  #   X0,
  #   ManifoldDiff.default_differential_backend()
  # )

## --------------------------------------------------------------------------------------------------------------
## ForwardDiff jacobian for Riemannian Levenberg-Marquardt
## --------------------------------------------------------------------------------------------------------------

struct JacF_RLM_ForwardDiff!{CF, JC}
  costF!::CF
  X0::Vector{Float64}
  res::Vector{Float64}
  sparsity::Union{SparseMatrixCSC, Nothing}
  jac_cache::JC  # ForwardColorJacCache or nothing
end

function JacF_RLM_ForwardDiff!(M, costF!, p, fg=nothing;
  all_points=p,
  basis_domain::AbstractBasis = LieGroups.DefaultLieAlgebraOrthogonalBasis(),
  is_sparse=!isnothing(fg),
)
  res = reduce(vcat, map(f -> f(all_points), Vector(costF!.costfuns)))
  X0 = zeros(manifold_dimension(M))

  if is_sparse && !isnothing(fg)
    factLabels = collect(getproperty.(costF!.costfuns, :faclbl))
    sparsity = eltype(res).(getSparsityPattern(fg, costF!.varLabels, factLabels))
    colorvec = matrix_colors(sparsity)
    # build the in-place wrapper that ForwardColorJacCache expects: f(out, x)
    function _inplace_costf!(out, Xc)
      _X = get_vector(M, p, Xc, basis_domain)
      _q = exp(M, p, _X)
      costF!(M, out, _q)
    end
    jac_cache = ForwardColorJacCache(_inplace_costf!, X0; dx=similar(res), colorvec, sparsity)
  else
    sparsity = nothing
    jac_cache = nothing
  end

  return JacF_RLM_ForwardDiff!(costF!, X0, res, sparsity, jac_cache)
end

function (jacF!::JacF_RLM_ForwardDiff!)(
  M::AbstractManifold,
  J,
  p;
  basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
  # basis_domain::AbstractBasis = LieGroups.DefaultLieAlgebraOrthogonalBasis(),
)
  X0 = jacF!.X0
  fill!(X0, 0)

  if !isnothing(jacF!.jac_cache)
    # sparse path: use coloring-aware ForwardDiff via SparseDiffTools
    function _inplace_costf_sparse!(out, Xc)
      X = get_vector(M, p, Xc, basis_domain)
      # X = hat(LieAlgebra(M), Xc)
      q = exp(M, p, X)
      jacF!.costF!(M, out, q)
    end
    forwarddiff_color_jacobian!(J, _inplace_costf_sparse!, X0, jacF!.jac_cache)
  else
    # dense path: standard ForwardDiff
    nres = length(jacF!.res)
    function costf(Xc)
      X = get_vector(M, p, Xc, basis_domain)
      q = exp(M, p, X)
      _res = zeros(eltype(Xc), nres)
      jacF!.costF!(M, _res, q)
      return _res
    end
    ForwardDiff.jacobian!(J, costf, X0)
  end
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

  # vec(CartesianIndices((R_range[2], C_range[1])))

  return sparse(getindex.(iter,1), getindex.(iter,2), ones(Bool, length(iter)))
end

function precisionFiniteDiff(M, jacF!::JacF_RLM!, p0)
    # Jcache
    X0 = fill!(deepcopy(jacF!.X0), 0)
    
    function costf(Xc)
      let res = jacF!.res, X = jacF!.X, q = jacF!.q, p0=p0
        get_vector!(M, X, p0, Xc, DefaultOrthogonalBasis())
        # get_vector!(M, X, p0, Xc, LieGroups.DefaultLieAlgebraOrthogonalBasis())
        exp!(M, q, p0, X)
        1/2*norm(jacF!.costF!(M, res, q))^2
      end
    end
    
    FiniteDiff.finite_difference_hessian(costf, X0)
end

function precisionFiniteDiff(M, jacF!::JacF_RLM_ForwardDiff!, p0)
    X0 = fill!(copy(jacF!.X0), 0)
    nres = length(jacF!.res)
    
    function costf(Xc)
      # X = get_vector(M, p0, Xc, DefaultOrthogonalBasis())
      X = get_vector(M, p0, Xc, LieGroups.DefaultLieAlgebraOrthogonalBasis())
      q = exp(M, p0, X)
      _res = zeros(nres)
      jacF!.costF!(M, _res, q)
      return 1/2*norm(_res)^2
    end
    
    FiniteDiff.finite_difference_hessian(costf, X0)
end

function qr_linear_subsolver!(sk, JJ, grad_f_c)
  sk .= qr(JJ) \ grad_f_c
  return sk
end

function pinv_subsolver!(sk, JJ, grad_f_c)
  sk .= pinv(JJ) * grad_f_c
  return sk
end

"""
    DebugTension(dof; io=stdout)

Manopt `DebugAction` that prints the tension (reduced chi-squared) at each iteration.
Tension = 2*cost / dof, where `dof = N - M` (residual dimension minus state dimension).

Usage in `solve_RLM`:
```julia
dof = num_components - manifold_dimension(M)
solve_RLM(fg; debug = [:Iteration, " | ", DebugTension(dof), "\n", 1])
```
"""
mutable struct DebugTension <: Manopt.DebugAction
  dof::Int
  io::IO
  format::String
end
DebugTension(dof::Int; io::IO=stdout, format="tension: %.2f") = DebugTension(dof, io, format)

function (d::DebugTension)(p::Manopt.AbstractManoptProblem, st::Manopt.AbstractManoptSolverState, k::Int)
  if d.dof <= 0
    s = NaN
  else
    cost = Manopt.get_cost(p, Manopt.get_iterate(st))
    s = sqrt(2 * cost / d.dof)
  end
  Printf.format(d.io, Printf.Format(d.format), s)
  return nothing
end

# replace :tension symbol with DebugTension(dof) in a debug vector
_inject_tension(debug, dof) = map(x -> x === :tension ? DebugTension(dof) : x, debug)

function solve_RLM(
  fg,
  varlabels = ls(fg),
  faclabels = lsf(fg);
  is_sparse = true,
  finiteDiffCovariance = false,
  jacobian_method::Symbol = :finitediff,
  solveKey::Symbol = :parametric,
  # linear_subsolver! = Manopt.default_lm_lin_solve!,
  linear_subsolver! = qr_linear_subsolver!,
  kwargs...
)

  # get the manifold and variable types
  vars = getVariable.(fg, varlabels)
   
  M, varTypes, vartypeslist = buildGraphSolveManifold(vars)

  varIntLabel, varlabelsAP = getVarIntLabelMap(vartypeslist)

  #Can use varIntLabel (because its an OrderedDict), but varLabelsAP makes the ArrayPartition.
  p0 = map(varlabelsAP) do label
    mean(getBelief(getState(fg, label, solveKey)))
  end

  # create an ArrayPartition{CalcFactorResidual} for faclabels
  calcfacs = CalcFactorResidualAP(fg, faclabels, varIntLabel)

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostFres!(calcfacs, collect(varlabelsAP))

  # jacobian of function for Riemannian Levenberg-Marquardt
  if jacobian_method == :forwarddiff
    jacF! = JacF_RLM_ForwardDiff!(M, costF!, p0, fg; is_sparse)
  else
    jacF! = JacF_RLM!(M, costF!, p0, fg; is_sparse)
  end

  num_components = length(jacF!.res)
  initial_residual_values = zeros(num_components)

  # initial_jacobian_f not type stable, but function barrier so should be ok.
  initial_jacobian_f = if jacF! isa JacF_RLM! && is_sparse
    jacF!.Jcache.sparsity
  elseif jacF! isa JacF_RLM_ForwardDiff! && !isnothing(jacF!.sparsity)
    jacF!.sparsity
  else
    zeros(num_components, manifold_dimension(M))
  end

  # inject DebugTension for :tension symbol in debug kwarg
  dof = num_components - manifold_dimension(M)
  if haskey(kwargs, :debug)
    kwargs = (; kwargs..., debug = _inject_tension(kwargs[:debug], dof))
  end

  lm_r = Manopt.LevenbergMarquardt!(
    M,
    costF!,
    jacF!,
    p0,
    num_components;
    evaluation=InplaceEvaluation(),
    jacobian_tangent_basis = LieGroups.DefaultLieAlgebraOrthogonalBasis(),
    # jacobian_tangent_basis = DefaultOrthogonalBasis(),
    initial_residual_values,
    initial_jacobian_f,
    linear_subsolver!,
    kwargs...
  )

  if finiteDiffCovariance
    Λ = precisionFiniteDiff(M, jacF!, lm_r)
  else
    J = initial_jacobian_f
    jacF!(M, J, lm_r) # recompute J at solution point
    Λ = Symmetric(J'J) # approx Hessian = precision matrix
  end

  # tension (reduced chi-squared): ||r||^2 / (N-M) = 2*cost / (N-M)
  dof = num_components - manifold_dimension(M)
  final_res = zeros(num_components)
  costF!(M, final_res, lm_r)
  tension = dof > 0 ? sum(abs2, final_res) / dof : NaN

  return M, varlabelsAP, lm_r, Λ, tension
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

function build_costF_jacF(
    fg,
    varlabels = ls(fg),
    faclabels = lsf(fg);
    is_sparse = false,
)
  
  # get the manifold and variable types
  vars = getVariable.(fg, varlabels)
    
  M, varTypes, vartypeslist = buildGraphSolveManifold(vars)

  varIntLabel, varlabelsAP = getVarIntLabelMap(vartypeslist)

  #Can use varIntLabel (because its an OrderedDict), but varLabelsAP makes the ArrayPartition.
  p0 = map(varlabelsAP) do label
    getVal(fg, label, solveKey = :parametric)[1]
  end

  # create an ArrayPartition{CalcFactorResidual} for faclabels
  calcfacs = CalcFactorResidualAP(fg, faclabels, varIntLabel)

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostFres!(calcfacs, collect(varlabelsAP))

  # jacobian of function for Riemannian Levenberg-Marquardt
  jacF! = JacF_RLM!(M, costF!, p0, fg; is_sparse)
  
  return M, costF!, jacF!, p0
end

function solve_RLM_conditional(
  fg,
  frontals::Vector{Symbol} = ls(fg),
  separators::Vector{Symbol} = setdiff(ls(fg), frontals);
  is_sparse=false,
  finiteDiffCovariance=true,
  jacobian_method::Symbol = :finitediff,
  solveKey::Symbol = :parametric,
  linear_subsolver! = qr_linear_subsolver!,
  kwargs...
)
  is_sparse && error("Sparse solve_RLM_conditional not supported yet")

  # get the subgraph formed by all frontals, separators and fully connected factors
  varlabels = union(frontals, separators)
  _, faclabels = listNeighborhood(fg, varlabels, 1)

  filter!(faclabels) do fl
    return issubset(getVariableOrder(fg, fl), varlabels)
  end

  @assert !isempty(faclabels) "Empty factor set for graph with variables $(ls(fg))"

  frontal_vars = getVariable.(fg, frontals)
  separator_vars = getVariable.(fg, separators)

  # so the subgraph consists of varlabels(frontals + separators) and faclabels

  _, _, frontal_vartypeslist = getVariableTypesCount(getVariable.(fg,frontals))
  frontal_varIntLabel, frontal_varlabelsAP = getVarIntLabelMap(frontal_vartypeslist)

  if isempty(separators)
    separator_vartypeslist = OrderedDict{DataType, Vector{Symbol}}()
    separator_varlabelsAP = ArrayPartition{Symbol,Tuple}(())
  else
    _, _, separator_vartypeslist = getVariableTypesCount(getVariable.(fg,separators))
    separator_varIntLabel, separator_varlabelsAP = getVarIntLabelMap(separator_vartypeslist)
  end

  all_varlabelsAP = ArrayPartition((frontal_varlabelsAP.x..., separator_varlabelsAP.x...))

  all_points = map(all_varlabelsAP) do label
    mean(getBelief(getState(fg, label, solveKey)))
  end
  
  p0 = ArrayPartition(all_points.x[1:length(frontal_varlabelsAP.x)])

  all_varIntLabel = OrderedDict{Symbol,Int}(
    map(enumerate(all_varlabelsAP)) do (i,l)
      l=>i
    end
  )
  # varIntLabel_frontals = filter(p->first(p) in frontals, varIntLabel)
  # varIntLabel_separators = filter(p->first(p) in separators, varIntLabel)

  calcfacs = CalcFactorResidualAP(fg, faclabels, all_varIntLabel)

  # get the manifold and variable types
   
  M, varTypes, vartypeslist = buildGraphSolveManifold(frontal_vars)
  
  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostFres_cond!(all_points, calcfacs, Vector{Symbol}(collect(all_varlabelsAP)))

  # jacobian of function for Riemannian Levenberg-Marquardt
  if jacobian_method == :forwarddiff
    jacF! = JacF_RLM_ForwardDiff!(M, costF!, p0, fg; all_points, is_sparse)
  else
    jacF! = JacF_RLM!(M, costF!, p0, fg; all_points, is_sparse)
  end

  num_components = length(jacF!.res)

  initial_residual_values = zeros(num_components)

  initial_jacobian_f = if jacF! isa JacF_RLM! && is_sparse
    jacF!.Jcache.sparsity
  elseif jacF! isa JacF_RLM_ForwardDiff! && !isnothing(jacF!.sparsity)
    jacF!.sparsity
  else
    zeros(num_components, manifold_dimension(M))
  end

  # inject DebugTension for :tension symbol in debug kwarg
  dof = num_components - manifold_dimension(M)
  if haskey(kwargs, :debug)
    kwargs = (; kwargs..., debug = _inject_tension(kwargs[:debug], dof))
  end

  lm_r = LevenbergMarquardt(
    M,
    costF!,
    jacF!,
    p0,
    num_components;
    evaluation=InplaceEvaluation(),
    initial_residual_values,
    initial_jacobian_f,
    linear_subsolver!,
    kwargs...
  )

  if finiteDiffCovariance
    Λ = precisionFiniteDiff(M, jacF!, lm_r)
  else
    jacF!(M, initial_jacobian_f, lm_r)
    Λ = Symmetric(initial_jacobian_f' * initial_jacobian_f)
  end

  # tension (reduced chi-squared): ||r||^2 / (N-M) = 2*cost / (N-M)
  final_res = zeros(num_components)
  costF!(M, final_res, lm_r)
  tension = sum(abs2, final_res) / dof
  
  return M, frontal_varlabelsAP, lm_r, Λ, tension
end

function extractMarginalsAP(M, labelsAP::ArrayPartition{Symbol}, Σ::AbstractArray{<:Real})
  st = 1
  Σvec = map(eachindex(labelsAP.x)) do i
      l = getDimension(M.manifolds[i].manifold)
      map(eachindex(labelsAP.x[i])) do j
          r = st:st + l - 1
          st += l
          SMatrix{l,l,Float64}(Σ[r,r])
      end
  end 
  ArrayPartition(Σvec...)
end

  #HEX solve
  # sparse J 0.025235 seconds (133.65 k allocations: 9.964 MiB
  # new1     0.013486 seconds (36.16 k allocations: 2.593 MiB)
  # new2    0.010764 seconds (34.61 k allocations: 3.111 MiB)
  # dense  J 0.022079 seconds (283.54 k allocations: 18.146 MiB)
  
function getInitOrderWavefront(fg, state_label::Symbol=:parametric; depth::Int=1)
    cliques = NamedTuple{(:frontals, :separators), Tuple{Vector{Symbol}, Vector{Symbol}}}[]
    
    all_vls = listVariables(fg)
    knowns = filter(vl -> hasState(fg, vl, state_label) && isInitialized(fg, vl, state_label), all_vls)
    unknowns = setdiff(all_vls, knowns)
    
    prior_vls, _ = listNeighborhood(fg, lsfPriors(fg), 1)
    
    while !isempty(unknowns)
        anchors = union(knowns, prior_vls)
        
        # Find Frontals
        neighborhood_vls, _ = isempty(anchors) ? (Symbol[], Symbol[]) : listNeighborhood(fg, anchors, 2 * depth)
        frontals = intersect(neighborhood_vls, unknowns)
        
        # Safety break if graph is completely floating (no priors, no knowns left to expand from)
        isempty(frontals) && break
        
        # Find Separators (Initialized variables touching our new frontals)
        frontal_neighbors, _ = listNeighborhood(fg, frontals, 2)
        separators = intersect(frontal_neighbors, knowns)
        
        # Record and Advance
        push!(cliques, (; frontals, separators))
        union!(knowns, frontals)
        setdiff!(unknowns, frontals)
    end
    
    return cliques
end

function autoinitParametric!(
  fg,
  clique_order = getInitOrderWavefront(fg; depth=3);
  reinit = false,
  kwargs...
)
  did_init = false
  @showprogress for cliq in clique_order
    did_init |= autoinitParametric!(fg, cliq.frontals, cliq.separators; reinit, kwargs...)
  end

  return did_init
end

function autoinitParametric!(dfg::AbstractDFG, initme::Symbol; kwargs...)
  return autoinitParametric!(dfg, getVariable(dfg, initme); kwargs...)
end

function autoinitParametric!(dfg::AbstractDFG, xi::VariableCompute; solveKey = :parametric, kwargs...)
  initme = getLabel(xi)
  prepareState!(xi, NLLSSolver(), solveKey)
  separators = ls2(dfg, initme)
  filter!(separators) do vl
    return hasState(dfg, vl, solveKey) && isInitialized(dfg, vl, solveKey)
  end
  return autoinitParametric!(dfg, [initme], separators; solveKey, kwargs...)
end

function autoinitParametric!(
  dfg::AbstractDFG,
  frontals::Vector{Symbol},
  separators::Vector{Symbol} = Symbol[];
  solveKey = :parametric,
  reinit::Bool = false,
  linear_subsolver! = pinv_subsolver!,
  kwargs...,
)
  #TODO prepare only the relevant states.
  prepareStates!(dfg, NLLSSolver(), solveKey)
  # Filter to only uninitialized variables (unless reinit)
  to_init = if reinit
    frontals
  else
    filter(v -> !isInitialized(dfg, v, solveKey), frontals)
  end
  isempty(to_init) && return false

  # Filter separators to only those already initialized
  active_separators = filter(separators) do vl
    hasState(dfg, vl, solveKey) && isInitialized(dfg, vl, solveKey)
  end

  # Nothing to initialize if no separators and no priors on any frontal
  if isempty(active_separators)
    has_any_prior = any(to_init) do v
      any(isPrior.(dfg, listNeighbors(dfg, v)))
    end
    has_any_prior || return false
  end

  # Check that we have usable factors
  varlabels = union(to_init, active_separators)
  _, faclabels = listNeighborhood(dfg, varlabels, 1)
  filter!(fl -> issubset(getVariableOrder(dfg, fl), varlabels), faclabels)
  isempty(faclabels) && return false

  # Prune stranded frontals
  connected_vars = unique(Iterators.flatten(getVariableOrder.(dfg, faclabels)))
  filter!(v -> v in connected_vars, to_init)
  isempty(to_init) && return false

  # Seed each frontal from an initialized separator of the same type
  for v in to_init
    xi = getVariable(dfg, v)
    vnd = getState(xi, solveKey)
    has_prior = any(isPrior.(dfg, listNeighbors(dfg, v)))
    if !has_prior && !isempty(active_separators)
      my_kind = getStateKind(xi)
      same_kind = filter(active_separators) do vl
        getStateKind(getVariable(dfg, vl)) === my_kind
      end
      if !isempty(same_kind)
        mn = mean(getBelief(getState(dfg, same_kind[1], solveKey)))
        _bw = cov(getBelief(vnd))
        # BW = getBW(getBelief(vnd))
        # _bw = (0<length(BW)) && isassigned(BW,1) ? BW[1] : nothing
        _hode = HomotopyDensity_legacy(getStateKind(vnd),[mn,]; bw=_bw, newbw=false)
        setBelief!(vnd, _hode)
      end
    end
    
    # perturb point slightly
    _M = getManifold(xi)
    tangent_coords = randn(manifold_dimension(_M)) * 1e-3
    X = get_vector(LieAlgebra(_M), tangent_coords)
    mn = mean(getBelief(vnd))
    mn_ = exp(_M, mn, X)
    bw = cov(getBelief(vnd))
    hode = HomotopyDensity_legacy(getStateKind(vnd), [mn_,]; bw, newbw=false)
    setBelief!(vnd, hode)
  end

  # Solve
  M, varlabelsAP, lm_r, Λ, _ = solve_RLM_conditional(dfg, to_init, active_separators; solveKey, linear_subsolver!, kwargs...)

  _Σ = (I)(size(Λ, 1))
  _Σ_ = sparse(1.0*I, size(Λ, 1), size(Λ, 1)) # create a sparse identity matrix
  # _Σ_ = (1.0*I)(size(Λ, 1)) # legacy was Bool, weird refactor forced premature Float64

  offset = 0

  invertonce = true
  # Update each frontal variable with result
  for (i, v) in enumerate(varlabelsAP)
    vrb = getVariable(dfg, v)
    state = getState(vrb, solveKey)

    # Update covariances from joint precision if positive definite
    if invertonce && !isnothing(Λ)
      F = cholesky!(Λ; check = false)
      if issuccess(F)
        invertonce = false
        _Σ_ .= (F \ _Σ) 
      end
    end
    dim = getDimension(vrb)
    r = (offset + 1):(offset + dim)
    offset += dim
    bw = ApproxManifoldProducts._forcestatic(_Σ_[r,r])
    # @info "WHAT" string(r) string(lm_r[i]) string(bw)
    hode = HomotopyDensity_legacy(getStateKind(state),[lm_r[i],]; bw, newbw=false)
    setBelief!(state, hode, true)
  end

  return true
end


"""
    $SIGNATURES

Batch parametric graph solve using Riemannian Levenberg Marquardt.
"""
solveGraphParametric(args...; kwargs...) = solve_RLM(args...; kwargs...)

function DFG.solveGraphParametric!(
  fg::AbstractDFG,
  args...; 
  init::Bool = true, 
  solveKey::Symbol = :parametric,
  is_sparse = true,
  # debug, stopping_criterion, damping_term_min=1e-2, 
  # expect_zero_residual=true,
  kwargs...
)
  # make sure variables has solverData, see #1637
  prepare!(fg, NLLSSolver(), solveKey)
  init && autoinitParametric!(fg; solveKey)

  M, v, r, Λ, tension = solve_RLM(fg, args...; is_sparse, kwargs...)

  updateParametricSolution!(fg, M, v, r, Λ)

  return M, v, r, Λ
end


## Check when time and delete if it can't be improved, curretnly ArrayPartition works best
#=
using FunctionWrappers: FunctionWrapper

# call with 
calcfacs = CalcFactorResidualWrapper(fg, faclabels, varIntLabel, all_points)
costF! = CostF_RLM_WRAP2!(all_points, calcfacs, map(cfm->size(cfm.obj.x.iΣ,1), calcfacs))

function CalcFactorResidualWrapper(fg, factorLabels::Vector{Symbol}, varIntLabel::OrderedDict{Symbol, Int64}, points::ArrayPartition)
  factypes, typedict, alltypes = getFactorTypesCount(getFactor.(fg, factorLabels))
  
  # skip non-numeric prior (MetaPrior)
  #TODO test... remove MetaPrior{T} something like this
  metaPriorKeys = filter(k->contains(string(k), "MetaPrior"), collect(keys(alltypes)))
  delete!.(Ref(alltypes), metaPriorKeys)

  calcfacs = map(factorLabels) do labels
    fct = getFactor(fg, labels)
    # moet ek 'n view in p0 in maak wat jy net p0 update en CFM view automaties, toets dit...
    cfm = IIF.CalcFactorResidual(fg, fct, varIntLabel, points)
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
