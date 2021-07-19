
# PoC on jacobians for a factor

using IncrementalInference
using TensorCast


##


# T_pt_args[:] = [(T1::Type{<:InferenceVariable}, point1); ...]
function _prepFactorGradientLambdas(fct::Union{<:AbstractRelativeMinimize,<:AbstractRelativeRoots}, 
                                    measurement::Tuple,
                                    T_pt_args...;
                                    h::Real=1e-8  ) # numerical diff perturbation size
  #
  # get manifolds for all variables
  M = (x->getManifold(x[1])).(T_pt_args)
  # use the temporary factor graph throughout
  # FORCED TO START AT EITHER :x0 or :x1
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
  f_dsi_h =       (d,s,i) -> IIF._evalFactorTemporary!(fct, d, measurement, arr_T_pth_s_i(s,i)..., tfg=tfg, newFactor=false, currNumber=0, _slack=slack_resid )
  # standard calculus derivative definition (in coordinate space)
  Δf_dsi =        (d,s,i, crd=coord_(d)) -> (@info "HERE" string(crd) string(f_dsi_h(d,s,i)) ; (f_dsi_h(d,s,i) - crd)./h)
  # jacobian block per s, for each i
  ▽f_ds =         (d,s, crd=coord_(d)) -> (@show crd; ((i)->Δf_dsi(d,s,i,crd)).(1:length(crd)))
  # jacobian stored in user provided matrix
  ▽f_ds_J! =      (J::AbstractMatrix{<:Real}, d, s) -> (J_ = ▽f_ds(d,s); (@cast J[_d,_s] = J_[_s][_d]))

  # TODO generalize beyond binary
  λ_fncs = () # Vector{Vector{Function}}()    # by factor's variable order
  
  # number of blocks
  nblks = length(T_pt_args)
  # length of all coordinate dimensions together
  λ_sizes = length.(coord_.(1:nblks))          # by factor's variable order
  len = sum(λ_sizes)
  # full jacobian matrix
  J_f = zeros(len, len)

  # build final lambdas which are mapped to the blocks of the full jacobian gradients matrix J_f
  Σd = 0 
  # each variable T, go down to coords, add eps to a coord, back to point and look at the change in residual (assumed in coords for AbstractRelative[Minimize/Roots])
  # TODO change `a_` to `s_` as variable selection by factor order
  for (d_, Tpt_d) in enumerate(T_pt_args)
    λ_row = ()
    len_d =λ_sizes[d_]
    Σs = 0
    for (s_, Tpt_s) in enumerate(T_pt_args)
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
  return J_f, λ_fncs, λ_sizes
end


##

# function _makeTupleGradientLambdas(λ_fncs, λ_sizes)
#   # number of blocks
#   nblks = length(λ_sizes)
#   # length of all coordinate dimensions together
#   len = sum(λ_sizes)
#   # full jacobian matrix
#   J_f = zeros(len, len)

#   # prep a new nested-array to convert into nested-tuple later 
#   totuple = Vector{Vector{Function}}(undef, nblks)
#   for i in 1:length(totuple)
#     totuple[i] = Vector{Function}(undef, nblks)
#   end
  
#   # build final lambdas which are mapped to the blocks of the full jacobian gradients matrix J_f
#   Σd = 0 
#   for (d_,len_d) in enumerate(λ_sizes)
#     Σs = 0
#     for (s_,len_s) in enumerate(λ_sizes)
#       # create a view into the full jacobian matrix at (d,s)
#       _J_ds = view(J_f, (1+Σd):(Σd+len_d), (1+Σs):(Σs+len_s))
#       totuple[d_][s_] = () -> (λ_fncs[d_][s_])(_J_ds)
#       Σs += len_s
#     end
#     Σd += len_d
#   end

#   # convert nested-vector to a nested-tuple for better type stability when used after this function
#   firsttuples = []
#   for tt in totuple
#     push!(firsttuples, tuple(tt...))
#   end
#   nestedtuple = tuple(firsttuples...)

#   # return the nested-tuple-lambdas and gradients matrix location
#   return nestedtuple, J_f, totuple
# end



##

pp = LinearRelative(MvNormal([10;0],[1 0; 0 1]))
measurement = ([10.0;0.0],)
T_pt_vec = [(ContinuousEuclid{2},[0;0.0]); (ContinuousEuclid{2},[9.5;0])]

##

J__, λ_fncs, λ_sizes = _prepFactorGradientLambdas(pp, measurement, T_pt_vec...);

##

# nt, J_f, tt = _makeTupleGradientLambdas(λ_fncs, λ_sizes);

##



##





#