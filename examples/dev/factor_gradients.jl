
# PoC on jacobians for a factor

using IncrementalInference

##

pp = LinearRelative(MvNormal([10;0],[1 0; 0 1]))


##


testFactorResidualBinary(pp, ([1;0],), (ContinuousEuclid{2}, [0;0]), (ContinuousEuclid{2},[0;0]) )


# T_pt_args[:] = [(T1::Type{<:InferenceVariable}, point1); ...]
function _prepFactorGradientLambdas(fct::Union{<:AbstractRelativeMinimize,<:AbstractRelativeRoots}, 
                                    measurement::Tuple,
                                    T_pt_args...;
                                    h::Real=1e-8  ) # numerical diff perturbation size
  #
  # get manifolds for all variables
  M = (x->getManifold(x[1])).(T_pt_args)
  # use the temporary factor graph throughout
  tfg = initfg()

  # gradients relative to coords requires 
  slack_resid = calcFactorResidualTemporary(fct, measurement, T_pt_args..., tfg=tfg)
  # coord_s = IIF._evalFactorTemporary!(fct, sfidx, measurement, T_pt_args..., tfg=tfg, newFactor=false, _slack=slack_resid )
  
  # perturb the coords of one variable on the factor, and return a new vector of TypePoint tuples
  # TODO, replace with retract operations instead
  coord_ =        (s) -> AMP.makeCoordsFromPoint(M[s], T_pt_args[s][2])
  coord_h =       (s,i, crd=coord_(s)) -> (crd[i] += h; crd)
  # reassemble TypePoint vector with perturbation at (s,i)
  T_pth_s_i =     (s,i) -> (T_pt_args[s][1], AMP.makePointFromCoords(M[s], coord_h(s,i), T_pt_args[s][2]) )
  arr_T_pth_s_i = (s,i) -> ((T_pt_args[1:(s-1)])..., T_pth_s_i(s,i), (T_pt_args[(s+1):end])...)
  # build a residual calculation specifically considering graph factor selections `s`, e.g. for binary `s ∈ {1,2}``.
  f_sih =         (s,i) -> calcFactorResidualTemporary(fct, measurement, arr_T_pth_s_i(s,i)..., tfg=tfg, newFactor=false, _slack=slack_resid )
  # standard calculus derivative definition (in coordinate space)
  Δf_si =         (s,i, crd=coord_(s)) -> (f_sih(s,i) - crd)./h

  # TODO generalize beyond binary
  λ_fncs = []    # by factor's variable order
  λ_sizes = []   # by factor's variable order

  # each variable T, go down to coords, add eps to a coord, back to point and look at the change in residual (assumed in coords for AbstractRelative[Minimize/Roots])
  # TODO change `a_` to `s_` as variable selection by factor order
  for (s_, Tpt_) in enumerate(T_pt_args)
    # assume binary factor to start, these are variables A and B
    Ta = Tpt_[1]
    pt_s = Tpt_[2]
      
    crd_s = coord_(s_)
    len_s = length(crd_s)
    push!(λ_sizes, len_s)

    # jacobian block per s, for each i
    ▽f_s = (s, crd=coord_(s)) -> ((i->Δf_si(s,i,crd)).(1:length(crd)))
    # jacobian stored in user provided matrix
    ▽f_s_J! = (J::AbstractMatrix{<:Real}, s) -> (J_ = ▽f_s(s); (@cast J[i,s] = J_[j][s]))
    # function is ready but calculation for actual jacobian values must still be done
    push!(λ_fncs, ▽f_J!)
  end
  
  # build Tuple{Tuple{<:Function,size}...}
  
end

#





#