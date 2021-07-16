
# PoC on jacobians for a factor

using IncrementalInference

##

pp = LinearRelative(MvNormal([10;0],[1 0; 0 1]))


##


testFactorResidualBinary(pp, ([1;0],), (ContinuousEuclid{2}, [0;0]), (ContinuousEuclid{2},[0;0]) )


# T_pt_args[:] = [(T1::Type{<:InferenceVariable}, point1); ...]
function _calcGradientFactorBinary( fct::Union{<:AbstractRelativeMinimize,<:AbstractRelativeRoots}, 
                                    meas::Tuple,
                                    T_pt_args...;
                                    h::Real=1e-8 ) # numerical diff perturbation size
  #
  # gradients relative to coords requires 
  resid_ = testFactorResidualBinary(fct, meas, T_pt_args...)
  
  # build a residual calculation specifically considering graph factor selections `s`, e.g. for binary `s ∈ {1,2}``.
  δ_s = (s,T_pt) -> testFactorResidualBinary(fct, meas, (T_pt_args[1:(s-1)])..., T_pt, (T_pt_args[(s+1):end])...)
  # TODO generalize beyond binary

  # each variable T, go down to coords, add eps to a coord, back to point and look at the change in residual (assumed in coords for AbstractRelative[Minimize/Roots])
  for (a_, Tpt_) in enumerate(T_pt_args)
    # assume binary factor to start, these are variables A and B
    Ta = Tpt_[1]
    Ma = getManifold(Ta)
    pt_a = Tpt_[2]
    u_a = pt_a
      
    # replace with retract operations only
    coords_a = AMP.makeCoordsFromPoint(Ma, pt_a)
    # perturb coords
    pr_coord_ = (i) -> (ca_ = deepcopy(coords_a); ca_[i] += h; ca_)  
    pt_a_i = (i) -> AMP.makePointFromCoords(Ma,pr_coord_(i),u_a)
    resid_a_i = (s,i) -> δ_s(s, (Ta, pt_a_i(i)))
    # standard calculus derivative definition (in coordinate space)
    Δf_h = (i) -> (resid_a_i(a_,i) - resid_)./h
    # jacobian block in reverse indices
    ▽f = () -> ((i->Δf_h(i)).(1:length(coords_a)))
    # jacobian stored in user provided matrix
    ▽f_ij! = (J) -> (J_ = ▽f(); (@cast J[i,j] = J_[j][i]))

    # function is ready but calculation must still be done
  end
  


end



#