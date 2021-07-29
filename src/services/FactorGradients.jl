# utilities for calculating the gradient over factors

export getCoordSizes


# T_pt_args[:] = [(T1::Type{<:InferenceVariable}, point1); ...]
# FORCED TO START AT EITHER :x1
function _prepFactorGradientLambdas(fct::Union{<:AbstractRelativeMinimize,<:AbstractRelativeRoots}, 
                                    measurement::Tuple,
                                    varTypes::Tuple,
                                    pts::Tuple;
                                    tfg::AbstractDFG = initfg(),
                                    # gradients relative to coords requires 
                                    slack_resid = calcFactorResidualTemporary(fct, measurement, varTypes, pts, tfg=tfg),
                                    # numerical diff perturbation size
                                    h::Real=1e-4  ) 
  #
  # get manifolds for all variables
  M = getManifold.(varTypes)
  # use the temporary factor graph throughout

  # TODO, replace with retract operations instead -- i.e. closer to Manifolds.jl tangent representations
  coord_ =        (s) -> AMP.makeCoordsFromPoint(M[s], pts[s])                         # vee(M,..., log(M,...) ) 
  # perturb the coords of one variable on the factor
  coord_h =       (s,i, crd=coord_(s)) -> (crd[i] += h; crd)
  # reassemble TypePoint vector with perturbation at (s,i)
  T_pth_s_i =     (s,i) -> AMP.makePointFromCoords(M[s], coord_h(s,i), pts[s])         # exp(M,..., hat(M,...) )
  tup_pt_s_i_h =  (s,i) -> tuple(pts[1:(s-1)]..., T_pth_s_i(s,i), pts[(s+1):end]...)
  # build a residual calculation specifically considering graph factor selections `s`, e.g. for binary `s ∈ {1,2}`.
  f_dsi_h =       (d,s,i) -> IIF._evalFactorTemporary!(fct, d, measurement, varTypes, tup_pt_s_i_h(s,i); tfg=tfg, newFactor=false, currNumber=0, _slack=slack_resid )
  # standard calculus derivative definition (in coordinate space)
  Δf_dsi =        (d,s,i, crd=coord_(d)) -> (f_dsi_h(d,s,i)[1] - crd)./h
  # jacobian block per s, for each i
  ▽f_ds =         (d,s, crd=coord_(d)) -> ((i)->Δf_dsi(d,s,i,crd)).(1:length(crd))
  # jacobian stored in user provided matrix
  ▽f_ds_J! =      (J::AbstractMatrix{<:Real}, d, s) -> (J_ = ▽f_ds(d,s); (@cast J[_d,_s] = J_[_s][_d]))

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
    len_d =λ_sizes[d_]
    Σs = 0
    for (s_, T_s) in enumerate(varTypes)
      len_s = λ_sizes[s_]
      # create a view into the full jacobian matrix at (d,s)
      _J_ds = view(J_f, (1+Σd):(Σd+len_d), (1+Σs):(Σs+len_s))
      # function is ready for calculation but actual jacobian values must still be done
      λ_row = tuple(λ_row..., () -> ▽f_ds_J!(_J_ds, d_,s_))
      Σs += len_s
    end
    λ_fncs = tuple(λ_fncs..., λ_row)
    Σd += len_d
  end
  
  # full gradients jacobian matrix, nested-tuple of lambdas to update, and sizes of blocks
  return J_f, λ_fncs, λ_sizes, slack_resid
end


function FactorGradientsCached!(fct::Union{<:AbstractRelativeMinimize, <:AbstractRelativeRoots},
                                varTypes::Tuple,
                                meas_single::Tuple, 
                                pts::Tuple; 
                                h::Real=1e-4  )
  #
  # working memory location for computations
  tfg = initfg()
  
  # permanent location for points and later reference
  # generate the necessary lambdas
  J__, λ_fncs, λ_sizes, slack_resid = _prepFactorGradientLambdas(fct, meas_single, varTypes, pts; tfg=tfg, h=1e-4);

  # get the one factor in tfg
  fctsyms = lsf(tfg)
  @assert length(fctsyms) == 1 "Expecting only a single factor in tfg"

  # generate an object containing all the machinery necessary more rapid factor gradients, see DevNotes for future improvements
  FactorGradientsCached!( getFactor(tfg, fctsyms[1]),
                          J__,
                          slack_resid,
                          meas_single,
                          pts,
                          tfg,
                          λ_fncs,
                          λ_sizes )
end


# Gradient matrix has individual blocks
getCoordSizes(fgc::FactorGradientsCached!) = fgc._coord_sizes


function (fgc::FactorGradientsCached!)(meas_pts...)
  # separate the measurements (forst) from the variable points (rest)
  lenm = length(fgc.measurement)
  @assert (length(fgc.currentPoints)+lenm) == length(meas_pts) "Unexpected number of arguments, got $(length(meas_pts)) but expected $(length(fgc.currentPoints)+lenm) instead.  Retry call with args (meas..., pts...)"

  # update in-place the new measurement value in preparation for new gradient calculation
  for (m, tup_m) in enumerate(fgc.measurement)
    _setPointsMani!(tup_m, meas_pts[m])
  end

  # update the residual _slack in preparation for new gradient calculation
  fct = getFactorType(fgc.dfgfct)
  measurement = tuple(meas_pts[1:lenm]...)
  pts = tuple(meas_pts[(1+lenm):end]...)
  varTypes = tuple(getVariableType.(getVariable.(fgc._tfg, getVariableOrder(fgc.dfgfct)))...)
  new_slack = calcFactorResidualTemporary(fct, measurement, varTypes, pts; tfg=fgc._tfg)
  # TODO make sure slack_residual is properly wired up with all the lambda functions as expected
  _setPointsMani!(fgc.slack_residual, new_slack)

  # set new points in preparation for new gradient calculation
  for (s,pt) in enumerate(meas_pts[(lenm+1):end])
    # update the local memory in fgc to take the values of incoming `pts`
    _setPointsMani!(fgc.currentPoints[s], pt)
  end

  # update the gradients at new values contained in fgc
  st = 0
  for (s,λ_tup) in enumerate(fgc._λ_fncs), (k,λ) in enumerate(λ_tup)
    # updating of the cached gradients assume that diagonal is zero for eigen value=1 -- i.e. (dA/dA - I)=0
    if s == k
      # coord length/size of this particular block
      szi = fgc._coord_sizes[s]
      _blk = view( fgc.cached_gradients, (st+1):(st+szi), (st+1):(st+szi) )
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