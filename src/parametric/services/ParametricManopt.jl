using Manopt
using FiniteDiff
using SparseDiffTools
using SparseArrays

# using ForwardDiff
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
  fct::DFGFactor, 
  varIntLabel
)
  fac_func = getFactorType(fct)
  varOrder = getVariableOrder(fct)

  varOrderIdxs = getindex.(Ref(varIntLabel), varOrder)

  M = getManifold(getFactorType(fct))

  dims = manifold_dimension(M)

  meas, iΣ = getFactorMeasurementParametric(fct)

  sqrt_iΣ = convert(SMatrix{dims, dims}, sqrt(iΣ))
  cache = preambleCache(fg, getVariable.(fg, varOrder), getFactorType(fct))

  return CalcFactorResidual(
    fct.label,
    getFactorMechanics(fac_func),
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
  p::AbstractArray{T},
  st::Int
) where {T, FT, N, D}
  for cfm in cfm_part
    x[st:st + D - 1] = cfm(p) #NOTE looks like do not broadcast here
    st += D
  end
  return st
end

function calcFactorResVec_threaded!(x::Vector{T}, cfm_part::Vector{<:CalcFactorResidual}, p::AbstractArray{T}, st::Int) where T
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

function (costf::CostFres!{CFT})(M::AbstractManifold, x::Vector{T}, p::AbstractVector{T}) where {CFT,T}
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
struct JacF_RLM!{CF, T, JC}
  costF!::CF
  X0::Vector{Float64}
  X::T
  q::T
  res::Vector{Float64}
  Jcache::JC
end

# function JacF_RLM!(M, costF!; basis_domain::AbstractBasis = DefaultOrthonormalBasis())
function JacF_RLM!(M, costF!, p, fg=nothing;
  all_points=p,
  basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
  is_sparse=!isnothing(fg)
)

  res = reduce(vcat, map(f -> f(all_points), Vector(costF!.costfuns)))

  X0 = zeros(manifold_dimension(M))
  
  X = get_vector(M, p, X0, basis_domain)

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
  basis_domain::AbstractBasis = DefaultOrthogonalBasis(),
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

# TODO only calculate marginal covariances

function covarianceFiniteDiff(M, jacF!::JacF_RLM!, p0)
    # Jcache
    X0 = fill!(deepcopy(jacF!.X0), 0)
    
    function costf(Xc)
      let res = jacF!.res, X = jacF!.X, q = jacF!.q, p0=p0
        get_vector!(M, X, p0, Xc, DefaultOrthogonalBasis())
        exp!(M, q, p0, X)
        1/2*norm(jacF!.costF!(M, res, q))^2
      end
    end
    
    H = FiniteDiff.finite_difference_hessian(costf, X0)

    # inv(H)
    Σ = try 
        Matrix(H) \ Matrix{eltype(H)}(I, size(H)...)
      catch ex #TODO only catch correct exception and try with pinv as fallback in certain cases.
        @warn "Hessian inverse failed" ex
        # Σ = pinv(H)
        nothing
      end
    return Σ
end

function solve_RLM(
  fg,
  varlabels = ls(fg),
  faclabels = lsf(fg);
  is_sparse = true,
  finiteDiffCovariance = false,
  solveKey::Symbol = :parametric,
  kwargs...
)

  # get the manifold and variable types
  vars = getVariable.(fg, varlabels)
   
  M, varTypes, vartypeslist = buildGraphSolveManifold(vars)

  varIntLabel, varlabelsAP = getVarIntLabelMap(vartypeslist)

  #Can use varIntLabel (because its an OrderedDict), but varLabelsAP makes the ArrayPartition.
  p0 = map(varlabelsAP) do label
    getVal(fg, label; solveKey)[1]
  end

  # create an ArrayPartition{CalcFactorResidual} for faclabels
  calcfacs = CalcFactorResidualAP(fg, faclabels, varIntLabel)

  #cost and jacobian functions
  # cost function f: M->ℝᵈ for Riemannian Levenberg-Marquardt 
  costF! = CostFres!(calcfacs, collect(varlabelsAP))

  # jacobian of function for Riemannian Levenberg-Marquardt
  jacF! = JacF_RLM!(M, costF!, p0, fg; is_sparse)

  num_components = length(jacF!.res)
  initial_residual_values = zeros(num_components)

  # initial_jacobian_f not type stable, but function barrier so should be ok.
  initial_jacobian_f = is_sparse ? 
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

  if length(initial_residual_values) < 1000 
    if finiteDiffCovariance
      # TODO this seems to be correct, but way to slow
      Σ = covarianceFiniteDiff(M, jacF!, lm_r)
    else
      # TODO make sure J initial_jacobian_f is updated, otherwise recalc jacF!(M, J, lm_r) # lm_r === p0
      J = initial_jacobian_f
      H = J'J # approx
      Σ = H \ Matrix{eltype(H)}(I, size(H)...)
      # Σ = pinv(H)
    end
  else
    @warn "Not estimating a Dense Covariance $(size(initial_jacobian_f))"
    Σ = nothing  
  end

  return M, varlabelsAP, lm_r, Σ
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

function solve_RLM_conditional(
  fg,
  frontals::Vector{Symbol} = ls(fg),
  separators::Vector{Symbol} = setdiff(ls(fg), frontals);
  is_sparse=false,
  finiteDiffCovariance=true,
  solveKey::Symbol = :parametric,
  kwargs...
)
  is_sparse && error("Sparse solve_RLM_conditional not supported yet")

  # get the subgraph formed by all frontals, separators and fully connected factors
  varlabels = union(frontals, separators)
  faclabels = sortDFG(setdiff(getNeighborhood(fg, varlabels, 1), varlabels))

  filter!(faclabels) do fl
    return issubset(getVariableOrder(fg, fl), varlabels)
  end

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
    seperator_varIntLabel, separator_varlabelsAP = getVarIntLabelMap(separator_vartypeslist)
  end

  all_varlabelsAP = ArrayPartition((frontal_varlabelsAP.x..., separator_varlabelsAP.x...))

  all_points = map(all_varlabelsAP) do label
    getVal(fg, label; solveKey)[1]
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
  jacF! = JacF_RLM!(M, costF!, p0, fg; all_points, is_sparse)

  num_components = length(jacF!.res)

  initial_residual_values = zeros(num_components)

  initial_jacobian_f = is_sparse ? 
    jacF!.Jcache.sparsity : 
    zeros(num_components, manifold_dimension(M))

  lm_r = LevenbergMarquardt(
    M,
    costF!,
    jacF!,
    p0,
    num_components;
    evaluation=InplaceEvaluation(),
    initial_residual_values,
    initial_jacobian_f,
    kwargs...
  )

  if finiteDiffCovariance
    Σ = covarianceFiniteDiff(M, jacF!, lm_r)
  else
    J = initial_jacobian_f
    Σ = pinv(J'J)
  end

  return M, all_varlabelsAP, lm_r, Σ
end

  #HEX solve
  # sparse J 0.025235 seconds (133.65 k allocations: 9.964 MiB
  # new1     0.013486 seconds (36.16 k allocations: 2.593 MiB)
  # new2    0.010764 seconds (34.61 k allocations: 3.111 MiB)
  # dense  J 0.022079 seconds (283.54 k allocations: 18.146 MiB)
  
function autoinitParametric!(
  fg,
  varorderIds = getInitOrderParametric(fg);
  reinit = false,
  kwargs...
)
  init_labels = @showprogress map(varorderIds) do vIdx
    autoinitParametric!(fg, vIdx; reinit, kwargs...)
  end
  filter!(!isnothing, init_labels)
  return init_labels
end

function autoinitParametric!(dfg::AbstractDFG, initme::Symbol; kwargs...)
  return autoinitParametric!(dfg, getVariable(dfg, initme); kwargs...)
end

function autoinitParametric!(
  dfg::AbstractDFG,
  xi::DFGVariable;
  solveKey = :parametric,
  reinit::Bool = false,
  perturb_point::Bool=false,
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
    
    # nothing to initialize if no initialized neighbors or priors
    if isempty(initfrom) && !any(isPrior.(dfg, listNeighbors(dfg, initme)))
      return false
    end

    vnd::VariableNodeData = getSolverData(xi, solveKey)
    
    if perturb_point
      _M = getManifold(xi)
      p = vnd.val[1]
      vnd.val[1] = exp(
        _M,
        p, 
        get_vector(
          _M,
          p,
          randn(manifold_dimension(_M))*10^-6,
          DefaultOrthogonalBasis()
        )
      )
    end
    M, vartypeslist, lm_r, Σ = solve_RLM_conditional(dfg, [initme], initfrom; solveKey, kwargs...)
    
    val = lm_r[1]
    vnd.val[1] = val

    !isnothing(Σ) && (vnd.bw .= Σ)
  
    # updateSolverDataParametric!(vnd, val, Σ)

    vnd.initialized = true
    #fill in ppe as mean
    Xc::Vector{Float64} = collect(getCoordinates(getVariableType(xi), val))
    ppe = MeanMaxPPE(solveKey, Xc, Xc, Xc)
    getPPEDict(xi)[solveKey] = ppe

    result = true

  else
    result = false
  end

  return result#isInitialized(xi, solveKey)
end


"""
    $SIGNATURES

Batch parametric graph solve using Riemannian Levenberg Marquardt.
"""
solveGraphParametric(args...; kwargs...) = solve_RLM(args...; kwargs...)

function DFG.solveGraphParametric!(
  fg::AbstractDFG,
  args...; 
  init::Bool = false, 
  solveKey::Symbol = :parametric,
  is_sparse = true,
  # debug, stopping_criterion, damping_term_min=1e-2, 
  # expect_zero_residual=true,
  kwargs...
)
  # make sure variables has solverData, see #1637
  makeSolverData!(fg; solveKey)
  if !(:parametric in fg.solverParams.algorithms)
    addParametricSolver!(fg; init = init)
  elseif init
    error("TODO: not implemented")
  end

  M, v, r, Σ = solve_RLM(fg, args...; is_sparse, kwargs...)

  updateParametricSolution!(fg, M, v, r, Σ)

  return M, v, r, Σ 
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

