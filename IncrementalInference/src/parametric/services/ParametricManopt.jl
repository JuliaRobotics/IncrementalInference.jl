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

"""
    $SIGNATURES

Coordinate index ranges per variable label, following the coordinate ordering of the
product manifold `M` and matching the `varlabelsAP` returned by [`getVarIntLabelMap`](@ref).

The container-level counterpart of `DensityBundlePoint.ranges`, for the RLM path that has a manifold
and a label partition but no bundle.
"""
function getCoordranges(M::AbstractManifold, varlabelsAP::ArrayPartition)
  ranges = OrderedDict{Symbol, UnitRange{Int}}()
  st = 1
  for i in eachindex(varlabelsAP.x)
    l = manifold_dimension(M.manifolds[i]) ÷ length(varlabelsAP.x[i])
    for j in eachindex(varlabelsAP.x[i])
      ranges[varlabelsAP.x[i][j]] = st:(st + l - 1)
      st += l
    end
  end
  return ranges
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
  metaPriorKeys = filter(k->contains(string(nameof(k)), "MetaPrior"), keys(alltypes))
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
Ideal value is 1.0, with values >> 1.0 indicating underfitting and < 1.0 indicating overfitting.

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
  linear_subsolver! = Manopt.default_lm_lin_solve!,
  kwargs...
)

  # get the manifold and variable types
  vars = getVariable.(fg, varlabels)
   
  M, varTypes, vartypeslist = buildGraphSolveManifold(vars)

  varIntLabel, varlabelsAP = getVarIntLabelMap(vartypeslist)

  #Can use varIntLabel (because its an OrderedDict), but varLabelsAP makes the ArrayPartition.
  p0 = map(varlabelsAP) do label
    getSingleModePoint(fg, label, solveKey)
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
    solveKey::Symbol = :parametric,
    p0 = nothing,
    partition = nothing,
)

  # get the manifold and variable types.
  M, varIntLabel, varlabelsAP = if isnothing(partition)
    vars = getVariable.(fg, varlabels)
    M_, _, vartypeslist = buildGraphSolveManifold(vars)
    vil, ap = getVarIntLabelMap(vartypeslist)
    (M_, vil, ap)
  else
    M_, ap, vil = buildPartitionedSolveManifold(map(g -> getVariable.(fg, g), partition))
    (M_, vil, ap)
  end

  #Can use varIntLabel (because its an OrderedDict), but varLabelsAP makes the ArrayPartition.
  if isnothing(p0)
    p0 = _readBasePoint(fg, varlabelsAP, solveKey)
  end

  # create an ArrayPartition{CalcFactorResidual} for faclabels
  calcfacs = CalcFactorResidualAP(fg, faclabels, varIntLabel)

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostFres!(calcfacs, collect(varlabelsAP))

  # jacobian of function for Riemannian Levenberg-Marquardt
  jacF! = JacF_RLM!(M, costF!, p0, fg; is_sparse)
  
  return M, costF!, jacF!, p0, varlabelsAP
end

"""
    _inflateFactorWeight(cf, fg, all_points, separator_set, solveKey)

Return `cf` with its whitening matrix `sqrt_iΣ` replaced by `sqrt(inv(Σ_eff))`, where the frozen
separator uncertainty has been folded into the factor's effective measurement covariance

    Σ_eff = Σ_fc + Σ_k J_sk Σ_sk J_skᵀ

`J_sk = ∂r/∂(separator k)` is evaluated once at the operating point `all_points`, by finite differences
of the unwhitened residual; per `dev/factor_jacobians.md` §4 a single evaluation suffices at `r ≈ 0`.
Factors between frontals only are returned unchanged.
"""
function _inflateFactorWeight(cf::CalcFactorResidual, fg, all_points, separator_set, solveKey::Symbol)
  sep_positions = findall(vl -> vl in separator_set, cf.varOrder)
  isempty(sep_positions) && return cf # frontal-only factor: weight unchanged

  D = getDimension(cf)
  points = map(idx -> all_points[idx], cf.varOrderIdxs)
  Σeff = inv(Symmetric(Matrix(cf.sqrt_iΣ' * cf.sqrt_iΣ))) # Σ_fc = inv(iΣ_fc)
  for pos in sep_positions
    svar = getVariable(fg, cf.varOrder[pos])
    G = getManifold(typeof(getStateKind(svar)))
    if G isa LieGroups.ValidationLieGroup
      G = G.lie_group #strip away ValidationLieGroup
    end
    Σ_s = Matrix{Float64}(getSingleModeCovariance(getState(svar, solveKey)))
    s0 = points[pos]
    # J_s = ∂(unwhitened residual)/∂(separator tangent), once, at the operating point
    Js = FiniteDiff.finite_difference_jacobian(zeros(manifold_dimension(G))) do Xc
      s = exp(G, s0, get_vector(G, s0, Xc, DefaultOrthogonalBasis()))
      cf(cf.meas, Base.setindex(points, s, pos)...)
    end
    Σeff = Σeff + Js * Σ_s * Js'
  end
  sqrt_iΣ_eff = convert(SMatrix{D, D}, sqrt(inv(Symmetric(Σeff))))
  return CalcFactorResidual(
    cf.faclbl, cf.factor, cf.varOrder, cf.varOrderIdxs, cf.meas, sqrt_iΣ_eff, cf.cache,
  )
end

"""
    _inflateSeparatorWeights(calcfacs, fg, all_points, separator_set, solveKey)

Copy the `calcfacs` `ArrayPartition`, inflating the whitening matrix of every factor that touches a
(frozen) separator via [`_inflateFactorWeight`](@ref). Preserves the partition/type structure.
"""
function _inflateSeparatorWeights(calcfacs, fg, all_points, separator_set, solveKey::Symbol)
  inflate(cf) = _inflateFactorWeight(cf, fg, all_points, separator_set, solveKey)
  parts = map(part -> map(inflate, part), calcfacs.x)
  return ArrayPartition{CalcFactorResidual, typeof(parts)}(parts)
end

## ================================================================================================
## Frontal solves: solve only `frontals`, with `separators` held at their means.
##
## Three probabilistic questions over one core, kept as separate functions because they carry
## different correctness contracts:
##
##   solve_RLM_conditional  P(F | S = μ_S)                       separators exact/frozen
##   solve_RLM_marginal     P(F) = ∫ P(F | S=s) P(S=s) ds        needs the separators' joint Σ_S
##   solve_RLM_propagate    joint MAP over (F,S) given S's prior  init/propagation; MOVES the mean
## ================================================================================================

"""
    $SIGNATURES

Internal shared core of [`solve_RLM_conditional`](@ref), [`solve_RLM_marginal`](@ref) and
[`solve_RLM_propagate`](@ref).  Solves `frontals` only, with `separators` frozen at their `solveKey`
means, and returns the solution plus the pieces each flavour needs to finish: the frozen operating
point `all_points`, the label partitions, the clique factor set, and the frontal-block precision
`Λ_FF`.

`reweight_separators=true` adds the second reweighted solve that makes an uncertain separator
*soften* rather than clamp its factors (see [`_inflateSeparatorWeights`](@ref)).  This changes the
point estimate and is only for [`solve_RLM_propagate`](@ref).
"""
function _solve_RLM_frontals_core(
  fg,
  frontals::Vector{Symbol},
  separators::Vector{Symbol};
  is_sparse=false,
  finiteDiffCovariance=true,
  jacobian_method::Symbol = :finitediff,
  solveKey::Symbol = :parametric,
  reweight_separators::Bool = false,
  linear_subsolver! = Manopt.default_lm_lin_solve!,
  kwargs...
)
  is_sparse && error("Sparse frontal solve not supported yet")

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
    getSingleModePoint(fg, label, solveKey)
  end

  p0 = ArrayPartition(all_points.x[1:length(frontal_varlabelsAP.x)])

  all_varIntLabel = OrderedDict{Symbol,Int}(
    map(enumerate(all_varlabelsAP)) do (i,l)
      l=>i
    end
  )

  calcfacs = CalcFactorResidualAP(fg, faclabels, all_varIntLabel)

  # get the manifold and variable types
  M, varTypes, vartypeslist = buildGraphSolveManifold(frontal_vars)

  # build cost + jacobian for a residual set and solve (separators stay frozen in `all_points`)
  function _build_and_solve(cfacs, p_start)
    costF! = CostFres_cond!(all_points, cfacs, Vector{Symbol}(collect(all_varlabelsAP)))
    if jacobian_method == :forwarddiff
      jacF! = JacF_RLM_ForwardDiff!(M, costF!, p_start, fg; all_points, is_sparse)
    else
      jacF! = JacF_RLM!(M, costF!, p_start, fg; all_points, is_sparse)
    end
    num_components = length(jacF!.res)
    initial_jacobian_f = if jacF! isa JacF_RLM! && is_sparse
      jacF!.Jcache.sparsity
    elseif jacF! isa JacF_RLM_ForwardDiff! && !isnothing(jacF!.sparsity)
      jacF!.sparsity
    else
      zeros(num_components, manifold_dimension(M))
    end
    dof = num_components - manifold_dimension(M)
    # inject DebugTension for :tension symbol in debug kwarg (local copy, safe to call twice)
    solve_kwargs = kwargs
    if haskey(solve_kwargs, :debug)
      solve_kwargs = (; solve_kwargs..., debug = _inject_tension(solve_kwargs[:debug], dof))
    end
    lm_r = LevenbergMarquardt(
      M,
      costF!,
      jacF!,
      p_start,
      num_components;
      evaluation=InplaceEvaluation(),
      initial_residual_values = zeros(num_components),
      initial_jacobian_f,
      linear_subsolver!,
      solve_kwargs...
    )
    return (; costF!, jacF!, lm_r, initial_jacobian_f, num_components, dof)
  end

  # 1. frozen solve: separators clamped at their means (the classic conditional solve)
  s = _build_and_solve(calcfacs, p0)

  # 2. optional reweight step (`solve_RLM_propagate`): at the operating point (r≈0) fold each frozen
  #    separator's covariance into its factor's effective noise, then 3. resolve reweighted so an
  #    uncertain separator *softens* (rather than clamps) its factor.  This yields the joint MAP over
  #    (F,S) and therefore MOVES the point estimate — off for the conditional/marginal inference path.
  #    See dev/factor_jacobians.md §4.
  if reweight_separators && !isempty(separators)
    all_points[1:length(s.lm_r)] .= s.lm_r
    calcfacs = _inflateSeparatorWeights(calcfacs, fg, all_points, Set(separators), solveKey)
    s = _build_and_solve(calcfacs, s.lm_r)
  end

  lm_r = s.lm_r

  # Frontal-block precision Λ_FF (separators frozen ⇒ this is the *conditional* precision).
  # With `reweight_separators`, the inflated `Σ_eff` weights already fold the separator uncertainty
  # into this ordinary frontal J'J — no separator block needed.
  if finiteDiffCovariance
    Λ_FF = precisionFiniteDiff(M, s.jacF!, lm_r)
  else
    s.jacF!(M, s.initial_jacobian_f, lm_r)
    Λ_FF = Symmetric(s.initial_jacobian_f' * s.initial_jacobian_f)
  end

  # tension (reduced chi-squared): ||r||^2 / (N-M) = 2*cost / (N-M)
  final_res = zeros(s.num_components)
  s.costF!(M, final_res, lm_r)
  tension = s.dof > 0 ? sum(abs2, final_res) / s.dof : NaN

  return (;
    M,
    frontal_varlabelsAP,
    lm_r,
    Λ_FF,
    tension,
    all_points,
    all_varlabelsAP,
    faclabels,
  )
end

"""
    $SIGNATURES

Conditional frontal solve: **`P(F | S = μ_S)`**, i.e. the separators are treated as *exactly* known
and clamped at their `solveKey` means.

Returns `(M, frontal_varlabelsAP, lm_r, Λ_FF, tension)` where `Λ_FF` is the **conditional** precision
of the frontals (it carries no separator uncertainty).

Use when the separators really are exact — e.g. the gauge anchor of a prior-free clique on the Bayes
tree upward pass, or a variable being initialised from already-fixed neighbours.  If the separators
carry uncertainty that should widen the frontals, use [`solve_RLM_marginal`](@ref) instead; if the
separators should be allowed to *move*, use [`solve_RLM_propagate`](@ref).
"""
function solve_RLM_conditional(
  fg,
  frontals::Vector{Symbol} = ls(fg),
  separators::Vector{Symbol} = setdiff(ls(fg), frontals);
  kwargs...
)
  c = _solve_RLM_frontals_core(fg, frontals, separators; kwargs...)
  return c.M, c.frontal_varlabelsAP, c.lm_r, c.Λ_FF, c.tension
end

"""
    $SIGNATURES

Marginal frontal solve: **`P(F) = ∫ P(F | S=s) P(S=s) ds`**, marginalising the separators out against
their supplied joint covariance `Σ_S`.

The point estimate is identical to [`solve_RLM_conditional`](@ref); only the covariance differs, by the
law of total covariance.  `Σ_S` must be the separators' **joint** covariance in the order of the
`separators` argument.  Returned `Σ_F`/`Σ_FS` follow the `frontals`/`separators` argument order, not
the type-grouped `frontal_varlabelsAP` order.

!!! note "Λ_SS is deliberately unused"
    Only the local `Λ_FF`/`Λ_FS` enter.  On the Bayes tree `Λ_SS` is the block that already left in
    this clique's upward message and is therefore already inside the parent's `Σ_S` — using it would
    double count.  That omission is what makes the downward pass correct with no cavity division.

Returns a NamedTuple `(; M, frontal_varlabelsAP, lm_r, Λ_FF, Λ_FS, Σ_F, Σ_FS, tension)`.
"""
function solve_RLM_marginal(
  fg,
  frontals::Vector{Symbol},
  separators::Vector{Symbol},
  Σ_S::AbstractMatrix;
  solveKey::Symbol = :parametric,
  finiteDiffCovariance = false,
  kwargs...
)
  c = _solve_RLM_frontals_core(
    fg, frontals, separators; solveKey, finiteDiffCovariance, kwargs...
  )

  # Joint [F;S] information at the solution.  Built over *all* clique variables so both blocks land
  # in one uniform tangent basis (the frontal-only jacF! of the core cannot supply Λ_FS).
  all_vars = Symbol[collect(c.all_varlabelsAP)...]
  _, _, all_vartypeslist = getVariableTypesCount(getVariable.(fg, all_vars))
  _, joint_varlabelsAP = getVarIntLabelMap(all_vartypeslist)

  # operating point: frontals at the solution, separators frozen at their means
  point_of = Dict{Symbol, Any}()
  for (i, lbl) in enumerate(collect(c.all_varlabelsAP))
    point_of[lbl] = c.all_points[i]
  end
  for i in eachindex(c.frontal_varlabelsAP.x), j in eachindex(c.frontal_varlabelsAP.x[i])
    point_of[c.frontal_varlabelsAP.x[i][j]] = c.lm_r.x[i][j]
  end

  M_all, _, jacF_all!, _ =
    build_costF_jacF(fg, all_vars, c.faclabels; is_sparse = false, solveKey)
  p_at = map(lbl -> point_of[lbl], joint_varlabelsAP)
  J = zeros(length(jacF_all!.res), manifold_dimension(M_all))
  jacF_all!(M_all, J, p_at)
  Λ_joint = J' * J

  ranges = getCoordranges(M_all, joint_varlabelsAP)
  idx_F = reduce(vcat, (collect(ranges[s]) for s in frontals))
  idx_S = reduce(vcat, (collect(ranges[s]) for s in separators))

  Λ_FF = Symmetric(Λ_joint[idx_F, idx_F])
  Λ_FS = Λ_joint[idx_F, idx_S]

  # marginalise S out — law of total covariance.  NOTE Λ_SS deliberately unused, see docstring.
  A = -(Λ_FF \ Λ_FS)
  Σ_F = Symmetric(A * Matrix(Σ_S) * A' + inv(Λ_FF))
  Σ_FS = A * Matrix(Σ_S)

  return (;
    c.M,
    c.frontal_varlabelsAP,
    c.lm_r,
    Λ_FF,
    Λ_FS,
    Σ_F,
    Σ_FS,
    c.tension,
  )
end

"""
    $SIGNATURES

Propagation / initialisation solve: the **joint MAP over `(F,S)`** given the separators' own
uncertainty, obtained by folding each frozen separator's covariance into its factors' effective noise
(`Σ_eff = Σ_fc + J_s Σ_s J_sᵀ`) and re-solving reweighted, so an uncertain separator *softens* rather
than clamps its factor.

Separator covariances are read per-variable from the graph's `solveKey` state, unlike
[`solve_RLM_marginal`](@ref), which takes an explicit joint `Σ_S`.

Returns `(M, frontal_varlabelsAP, lm_r, Λ, tension)`, with the separator uncertainty already in `Λ`.

!!! note "This moves the point estimate"
    The result is **not** the conditional mean and **not** the mean of `P(F)`.
"""
function solve_RLM_propagate(
  fg,
  frontals::Vector{Symbol} = ls(fg),
  separators::Vector{Symbol} = setdiff(ls(fg), frontals);
  kwargs...
)
  c = _solve_RLM_frontals_core(
    fg, frontals, separators; reweight_separators = true, kwargs...
  )
  return c.M, c.frontal_varlabelsAP, c.lm_r, c.Λ_FF, c.tension
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

"""
    getInitOrderSerial(fg, state_label=:parametric; depth=1)

Init order that advances **one variable at a time** — BFS from the prior-carrying variables outward,
each new frontal conditioned on its already-initialized neighbours.  Use via
`autoinitParametric!(fg, getInitOrderSerial(fg))`.

Prefer [`getInitOrderWavefront`](@ref) when variables are only *jointly* determined, e.g. partial priors,
where one variable cannot be initialized alone.
"""
function getInitOrderSerial(fg, state_label::Symbol=:parametric; depth::Int=1)
    cliques = NamedTuple{(:frontals, :separators), Tuple{Vector{Symbol}, Vector{Symbol}}}[]

    all_vls = listVariables(fg)
    knowns = filter(vl -> hasState(fg, vl, state_label) && isInitialized(fg, vl, state_label), all_vls)
    unknowns = setdiff(all_vls, knowns)

    # prior-carrying variables can be initialized on their own -> do them first, as anchors
    prior_neighbours, _ = listNeighborhood(fg, lsfPriors(fg), 1)
    prior_vls = intersect(prior_neighbours, unknowns)
    if !isempty(prior_vls)
        push!(cliques, (; frontals = prior_vls, separators = Symbol[]))
        union!(knowns, prior_vls)
        setdiff!(unknowns, prior_vls)
    end

    # then one variable per step, each conditioned on its already-initialized neighbours
    while !isempty(unknowns) && !isempty(knowns)
        ring_vls, _ = listNeighborhood(fg, knowns, 2 * depth)
        frontier = intersect(ring_vls, unknowns)
        isempty(frontier) && break

        v = first(frontier)
        v_neighbors, _ = listNeighborhood(fg, [v], 2)
        separators = intersect(v_neighbors, knowns)

        push!(cliques, (; frontals = [v], separators))
        push!(knowns, v)
        setdiff!(unknowns, [v])
    end

    return cliques
end

function autoinitParametric!(
  fg,
  clique_order = nothing;
  solveKey::Symbol = :parametric,
  reinit = false,
  kwargs...
)
  order = if isnothing(clique_order)
    getInitOrderWavefront(fg, solveKey; depth=3) #TODO maybe make default depth=1
  else
    clique_order
  end

  did_init = false
  @showprogress for cliq in order
    did_init |= autoinitParametric!(fg, cliq.frontals, cliq.separators; solveKey, reinit, kwargs...)
  end

  return did_init
end

function autoinitParametric!(dfg::AbstractDFG, initme::Symbol; kwargs...)
  return autoinitParametric!(dfg, getVariable(dfg, initme); kwargs...)
end

function autoinitParametric!(dfg::AbstractDFG, xi::VariableCompute; solveKey = :parametric, kwargs...)
  initme = getLabel(xi)
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
  linear_subsolver! = Manopt.default_lm_lin_solve!,
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
    # has_prior = any(isPrior.(dfg, listNeighbors(dfg, v)))
    # if !has_prior && !isempty(active_separators)
    if !isempty(active_separators)
      my_kind = getStateKind(xi)
      same_kind = filter(active_separators) do vl
        getStateKind(getVariable(dfg, vl)) === my_kind
      end
      if !isempty(same_kind)
        setSingleModeBelief!(vnd, getSingleModePoint(dfg, same_kind[1], solveKey); initialized = false)
      end
    end

    # perturb point slightly
    _M = getManifold(xi)
    tangent_coords = randn(manifold_dimension(_M)) * 1e-3
    mn = getSingleModePoint(vnd)
    X = hat(LieAlgebra(_M), tangent_coords, typeof(mn))
    setSingleModeBelief!(vnd, exp(_M, mn, X); initialized = false)
  end

  # Solve
  M, varlabelsAP, lm_r, Λ, _ = solve_RLM_conditional(dfg, to_init, active_separators; solveKey, linear_subsolver!, kwargs...)

  # Update each frontal variable with result
  for (i, v) in enumerate(varlabelsAP)
    vnd = getState(dfg, v, solveKey)
    setSingleModeBelief!(vnd, lm_r[i])
  end

  # Update covariances from joint precision if positive definite
  if !isnothing(Λ)
    F = cholesky!(Λ; check = false)
    if issuccess(F)
      Σ = F \ I(size(Λ, 1))
      offset = 0
      for (i, v) in enumerate(varlabelsAP)
        dim = manifold_dimension(getManifold(getVariable(dfg, v)))
        r = (offset + 1):(offset + dim)
        vnd = getState(dfg, v, solveKey)
        setSingleModeBelief!(vnd, getSingleModePoint(vnd), Σ[r, r])
        offset += dim
      end
    end
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

  M, v, r, Λ, tension = solve_RLM(fg, args...; solveKey, is_sparse, kwargs...)

  updateParametricSolution!(fg, M, v, r, Λ; solveKey)

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
  metaPriorKeys = filter(k->contains(string(nameof(k)), "MetaPrior"), keys(alltypes))
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
