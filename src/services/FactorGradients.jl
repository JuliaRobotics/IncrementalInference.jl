# utilities for calculating the gradient over factors



function factorJacobian(
  fg,
  faclabel::Symbol,
  p0 = ArrayPartition(first.(getVal.(fg, getVariableOrder(fg, faclabel), solveKey = :parametric))...),
  backend = ManifoldDiff.TangentDiffBackend(
    if v"0.4" <=  pkgversion(ManifoldDiff)
      ManifoldDiff.AutoFiniteDifferences(central_fdm(5, 1))
    else
      ManifoldDiff.FiniteDifferencesBackend()
    end
  )
)

  fac = getFactor(fg, faclabel)
  varlabels = getVariableOrder(fac)
  varIntLabel = OrderedDict(zip(varlabels, eachindex(varlabels)))

  cfm =  CalcFactorResidual(fg, fac, varIntLabel)
  
  function costf(p)
    points = map(idx->p.x[idx], cfm.varOrderIdxs)
    return cfm.sqrt_iΣ * cfm(cfm.meas, points...)
  end

  M_dom = ProductManifold(getManifold.(fg, varlabels)...)
  #TODO verify M_codom
  M_codom = Euclidean(manifold_dimension(getManifold(fac)))
  # Jx(M, p) = ManifoldDiff.jacobian(M, M_codom, calcfac, p, backend)
  
  return ManifoldDiff.jacobian(M_dom, M_codom, costf, p0, backend)
end

function factorGradient(
  cf::CalcFactorResidual,
  M,
  p,
  backend = ManifoldDiff.TangentDiffBackend(ManifoldDiff.FiniteDiffBackend()),
)
  ManifoldDiff.gradient(M, (x) -> 1//2 * norm(cf(x))^2, p, backend)
end

function factorJacobian(
  cf::CalcFactorResidual,
  M_dom,
  p,
  backend = ManifoldDiff.TangentDiffBackend(ManifoldDiff.FiniteDiffBackend()),
)
  # M_dom = ProductManifold(getManifold.(fg, varlabels)...)
  M_codom = Euclidean(manifold_dimension(getManifold(cf)))
    
  return ManifoldDiff.jacobian(M_dom, M_codom, cf, p, backend)
end

#
function factorGradient(
  cf::CalcFactorNormSq,
  M,
  p,
  backend = ManifoldDiff.TangentDiffBackend(ManifoldDiff.FiniteDiffBackend()),
)
  ManifoldDiff.gradient(M, cf, p, backend)
end

export getCoordSizes
export checkGradientsToleranceMask, calcPerturbationFromVariable

# T_pt_args[:] = [(T1::Type{<:InferenceVariable}, point1); ...]
# FORCED TO START AT EITHER :x1
function _prepFactorGradientLambdas(
  fct::Union{
    <:AbstractRelativeMinimize,
    # <:AbstractRelativeRoots,
    <:AbstractManifoldMinimize,
  },
  measurement,
  varTypes::Tuple,
  pts::Tuple;
  tfg::AbstractDFG = initfg(),
  _blockRecursion::Bool = true,
  # gradients relative to coords requires 
  slack_resid = calcFactorResidualTemporary(
    fct,
    varTypes,
    measurement,
    pts;
    tfg = tfg,
    _blockRecursion = _blockRecursion,
  ),
  # numerical diff perturbation size
  h::Real = 1e-4,
)
  #
  # get manifolds for all variables
  M = getManifold.(varTypes)
  # use the temporary factor graph throughout

  # TODO, replace with retract operations instead -- i.e. closer to Manifolds.jl tangent representations
  coord_ = (s) -> AMP.makeCoordsFromPoint(M[s], pts[s])                         # vee(M,..., log(M,...) ) 
  # perturb the coords of one variable on the factor
  coord_h = (s, i, crd = coord_(s)) -> (crd[i] += h; crd)
  # reassemble TypePoint vector with perturbation at (s,i)
  T_pth_s_i = (s, i) -> makePointFromCoords(M[s], coord_h(s, i), pts[s])         # exp(M,..., hat(M,...) )
  tup_pt_s_i_h = (s, i) -> tuple(pts[1:(s - 1)]..., T_pth_s_i(s, i), pts[(s + 1):end]...)
  # build a residual calculation specifically considering graph factor selections `s`, e.g. for binary `s ∈ {1,2}`.
  f_dsi_h =
    (d, s, i) -> IIF._evalFactorTemporary!(
      fct,
      varTypes,
      d,
      [measurement],
      tup_pt_s_i_h(s, i);
      tfg = tfg,
      newFactor = false,
      currNumber = 0,
      _slack = slack_resid,
    )
  # standard calculus derivative definition (in coordinate space)
  Δf_dsi =
    (d, s, i, crd = coord_(d)) ->
      (makeCoordsFromPoint(M[s], f_dsi_h(d, s, i)[1]) - crd) ./ h
  # jacobian block per s, for each i
  ▽f_ds = (d, s, crd = coord_(d)) -> ((i) -> Δf_dsi(d, s, i, crd)).(1:length(crd))
  # jacobian stored in user provided matrix
  ▽f_ds_J! =
    (J::AbstractMatrix{<:Real}, d, s) -> (J_ = ▽f_ds(d, s); (@cast J[_d, _s] = J_[_s][_d]))

  # TODO generalize beyond binary
  λ_fncs = ()                                  # by factor's variable order

  # number of blocks
  nblks = length(varTypes)
  # length of all coordinate dimensions together
  λ_sizes = length.(coord_.(1:nblks))          # by factor's variable order
  len = sum(λ_sizes)
  # full jacobian matrix
  J_f = zeros(len, len)

  # build final lambdas which are mapped to the blocks of the full jacobian gradients matrix J_f
  Σd = 0
  # each variable T, go down to coords, add eps to a coord, back to point and look at the change in residual (assumed in coords for AbstractRelative[Minimize/Roots])
  # TODO change `a_` to `s_` as variable selection by factor order
  for (d_, T_d) in enumerate(varTypes)
    λ_row = ()
    len_d = λ_sizes[d_]
    Σs = 0
    for (s_, T_s) in enumerate(varTypes)
      len_s = λ_sizes[s_]
      # create a view into the full jacobian matrix at (d,s)
      _J_ds = view(J_f, (1 + Σd):(Σd + len_d), (1 + Σs):(Σs + len_s))
      # function is ready for calculation but actual jacobian values must still be done
      λ_row = tuple(λ_row..., () -> ▽f_ds_J!(_J_ds, d_, s_))
      Σs += len_s
    end
    λ_fncs = tuple(λ_fncs..., λ_row)
    Σd += len_d
  end

  # full gradients jacobian matrix, nested-tuple of lambdas to update, and sizes of blocks
  return J_f, λ_fncs, λ_sizes, slack_resid
end

function FactorGradientsCached!(
  fct::Union{
    <:AbstractRelativeMinimize,
    <:AbstractManifoldMinimize,
  },
  varTypes::Tuple,
  meas_single,
  pts::Tuple;
  h::Real = 1e-4,
  _blockRecursion::Bool = true,
)
  #
  # working memory location for computations
  tfg = initfg()

  # permanent location for points and later reference
  # generate the necessary lambdas
  J__, λ_fncs, λ_sizes, slack_resid = _prepFactorGradientLambdas(
    fct,
    meas_single,
    varTypes,
    pts;
    tfg = tfg,
    _blockRecursion = _blockRecursion,
    h = 1e-4,
  )

  # get the one factor in tfg
  fctsyms = lsf(tfg)
  @assert length(fctsyms) == 1 "Expecting only a single factor in tfg"

  # generate an object containing all the machinery necessary more rapid factor gradients, see DevNotes for future improvements
  return FactorGradientsCached!(
    getFactor(tfg, fctsyms[1]),
    J__,
    slack_resid,
    meas_single,
    pts,
    tfg,
    λ_fncs,
    λ_sizes,
    h,
  )
end

# Gradient matrix has individual blocks
getCoordSizes(fgc::FactorGradientsCached!) = fgc._coord_sizes

function _setFGCSlack!(fgc::FactorGradientsCached!{F}, slack) where {F}
  return setPointsMani!(fgc.slack_residual, slack)
end

function _setFGCSlack!(
  fgc::FactorGradientsCached!{F, S},
  slack::Number,
) where {F, S <: Number}
  return fgc.slack_residual = slack
end

function (fgc::FactorGradientsCached!)(meas_pts...)

  # separate the measurements (forst) from the variable points (rest)
  lenm = 1# fgc.measurement is a single measurement so length is always 1

  @assert (length(fgc.currentPoints) + lenm) == length(meas_pts) "Unexpected number of arguments, got $(length(meas_pts)) but expected $(length(fgc.currentPoints)+lenm) instead.  Retry call with args (meas..., pts...), got meas_pts=$meas_pts"

  # update in-place the new measurement value in preparation for new gradient calculation
  # TODO should outside measurement be used or only that stored in FGC object?
  # for (m, tup_m) in enumerate(fgc.measurement)
  # setPointsMani!(tup_m, meas_pts[m])
  # end
  fgc.measurement = meas_pts[1] # why not 1:1 since ccwl.measurement::Vector{typeof(z)}

  # update the residual _slack in preparation for new gradient calculation
  fct = getFactorType(fgc.dfgfct)
  measurement = meas_pts[1]
  pts = meas_pts[2:end]
  varTypes =
    tuple(getVariableType.(getVariable.(fgc._tfg, getVariableOrder(fgc.dfgfct)))...)
  new_slack = calcFactorResidualTemporary(fct, varTypes, measurement, pts; tfg = fgc._tfg)
  # TODO make sure slack_residual is properly wired up with all the lambda functions as expected
  _setFGCSlack!(fgc, new_slack)
  # setPointsMani!(fgc.slack_residual, new_slack)

  # set new points in preparation for new gradient calculation
  for (s, pt) in enumerate(meas_pts[2:end])
    # update the local memory in fgc to take the values of incoming `pts`
    setPointsMani!(fgc.currentPoints[s], pt)
  end
  println.(fgc.currentPoints)

  # update the gradients at new values contained in fgc
  st = 0
  for (s, λ_tup) in enumerate(fgc._λ_fncs), (k, λ) in enumerate(λ_tup)
    # updating of the cached gradients assume that diagonal is zero for eigen value=1 -- i.e. (dA/dA - I)=0
    if s == k
      # coord length/size of this particular block
      szi = fgc._coord_sizes[s]
      _blk = view(fgc.cached_gradients, (st + 1):(st + szi), (st + 1):(st + szi))
      fill!(_blk, 0.0)
      # move on to next diagonal block
      st += szi
      continue
    end
    # recalculate the off diagonals
    λ()
  end

  # return newly calculated gradients
  return fgc.cached_gradients
end

# convenience function to update the gradients based on current measurement and point information stored in the fgc object
(fgc::FactorGradientsCached!)() = fgc(fgc.measurement, fgc.currentPoints...)

"""
    $SIGNATURES

Return a mask of same size as gradients matrix `J`, indicating which elements are above the expected sensitivity threshold `tol`.

Notes
- Threshold accuracy depends on two parts,
  - Numerical gradient perturbation size `fgc._h`,
  - Accuracy tolerance to which the factor residual is computed (not controlled here)
"""
function checkGradientsToleranceMask(
  fgc::FactorGradientsCached!,
  J::AbstractArray = fgc.cached_gradients;
  tol::Real = 0.02 * fgc._h,
)
  #
  # ignore anything 10 times smaller than numerical gradient delta used
  # NOTE this ignores the factor residual solve accuracy
  return tol * fgc._h .< abs.(J)
end

"""
    $SIGNATURES

Return a tuple of infoPerCoord vectors that result from input variables as vector of `::Pair`, i.e. `fromVar::Int => infoPerCoord`.
For example, a binary `LinearRelative` factor has a one-to-one influence from the input to the one other variable.

Notes

- Assumes the gradients in `fgc` are up to date -- if not, first run `fgc(measurement..., pts...)`.
- `tol` does not recalculate the gradients to a new tolerance, instead uses the cached value in `fgc` to predict accuracy.

Example

```julia
# setup
fct = LinearRelative(MvNormal([10;0],[1 0; 0 1]))
measurement = ([10.0;0.0],)
varTypes = (ContinuousEuclid{2}, ContinuousEuclid{2})
pts = ([0;0.0], [9.5;0])

# create the gradients functor object
fgc = FactorGradientsCached!(fct, varTypes, measurement, pts);
# must first update the cached gradients
fgc(measurement..., pts...)

# check the perturbation influence through gradients on factor
ret = calcPerturbationFromVariable(fgc, [1=>[1;1]])

@assert isapprox(ret[2], [1;1])
```

DevNotes
- FIXME Support n-ary source factors by extending `fromVar` to more than just one. 

Related

[`FactorGradientsCached!`](@ref), [`checkGradientsToleranceMask`](@ref)
"""
function calcPerturbationFromVariable(
  fgc::FactorGradientsCached!,
  from_var_ipc::AbstractVector{<:Pair};
  tol::Real = 0.02 * fgc._h,
)
  #
  blkszs = getCoordSizes(fgc)
  # assume projection through pp-factor from first to second variable 
  # ipc values from first variable belief, and zero for second
  ipcAll = zeros(sum(blkszs))

  # set any incoming infoPerCoord values
  for (fromVar, infoPC) in from_var_ipc
    # check on sizes with print warning
    if (blkszs[fromVar] == length(infoPC))
      nothing
    else
      @warn(
      "Expecting incoming length(infoPerCoord) to equal the block size for variable $fromVar, as per factor used to construct the FactorGradientsCached!: $(getFactorType(fgc.dfgfct))"
    )
    end
    # get range of interest
    curr_b = sum(blkszs[1:(fromVar - 1)]) + 1
    curr_e = sum(blkszs[1:fromVar])
    ipcAll[curr_b:curr_e] .= infoPC
  end

  # clamp gradients below numerical solver resolution
  mask = checkGradientsToleranceMask(fgc; tol = tol)
  J = fgc.cached_gradients
  _J = zeros(size(J)...)
  _J[mask] .= J[mask]

  # calculate the gradient influence on other variables
  ipc_pert = _J * ipcAll

  # round up over numerical solution tolerance
  dig = floor(Int, log10(1 / tol))
  ipc_pert .= round.(ipc_pert, digits = dig)

  # slice the result
  ipcBlk = []
  for (i, sz) in enumerate(blkszs)
    curr_b = sum(blkszs[1:(i - 1)]) + 1
    curr_e = sum(blkszs[1:i])
    blk_ = view(ipc_pert, curr_b:curr_e)
    push!(ipcBlk, blk_)
  end

  return tuple(ipcBlk...)
end

function calcPerturbationFromVariable(
  ccwl::CommonConvWrapper,
  sfidx::Int,
  smpid::Int = 1;
  tol::Real = 0.02 * ccwl.gradients_cached._h,
)
  #
  # collapse the hypo associated with smpid
  # get the variables associated with this hypo
  # assemble the leave-one-out of varidx=>infoPerCoords -- e.g. sfidx=1, `var_ipcs::Vector{<:Pair}=[2=>ipc2;]`
  #  NOTE varidx as per the factor args, i.e. after fractional associations (hypos) are collapsed
  # calcPerturbationFromVariable(ccwl.gradients_cached, var_ipcs; tol=tol)
  return error("UNDER CONSTRUCTION")
end

#
