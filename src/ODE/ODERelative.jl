
@info "IncrementalInference.jl is loading tools related to DifferentialEquations.jl"

using .DifferentialEquations

import IncrementalInference: getSample

export ODERelative


"""
    $TYPEDEF

Build a full ODE solution into a relative factor to condense possible sensor data into a relative transformation,
but keeping the parameter estimation process fluid.  Assumes first and second variable in order 
are of same dimension and compatible manifolds, such that ODE runs from Xi to Xi+1 on all
dimensions.  Internal state vector can be decoupled onto different domain as needed.

Notes
- Based on DifferentialEquations.jl
- `getSample` step does the `solve(ODEProblem)` step.
- `tspan` is taken from variables only once at object construction -- i.e. won't detect changed timestamps.
- Regular factor evaluation is done as full dimension `AbstractRelativeRoots`, and is basic linear difference.

DevNotes
- FIXME see 1025, `multihypo=` will not yet work. 
- FIXME Lots of consolidation and standardization to do, see RoME.jl #244 regarding Manifolds.jl.
- TODO does not yet handle case where a factor spans across two timezones.
"""
struct ODERelative{T <:InferenceVariable, P, D} <: AbstractRelativeRoots
  domain::Type{T}
  forwardProblem::P
  backwardProblem::P
  # second element of this data tuple is additional variables that will be passed down as a parameter
  data::D
  specialSampler::Function
end


"""
$SIGNATURES

Calculate a DifferentialEquations.jl ready `tspan::Tuple{Float64,Float64}` from DFGVariables.

DevNotes
- TODO does not yet incorporate Xi.nanosecond field.
- TODO does not handle timezone crossing properly yet.
"""
function _calcTimespan( Xi::AbstractVector{<:DFGVariable}  )
  #
  tsmps = getTimestamp.(Xi[1:2]) .|> DateTime .|> datetime2unix
  # toffs = (tsmps .- tsmps[1]) .|> x-> elemType(x.value*1e-3)
  return (tsmps...,)
end
# Notes
# - Can change numerical data return type using an additional first argument, `_calcTimespan(Float32, Xi)`.
# _calcTimespan(Xi::AbstractVector{<:DFGVariable}) = _calcTimespan(Float64, Xi)

# performance helper function
_maketuplebeyond2args = (w1=nothing,w2=nothing,w3_...) -> (w3_...,)


function ODERelative( Xi::AbstractVector{<:DFGVariable},
                      domain::Type{<:InferenceVariable},
                      f::Function,
                      data=()->();
                      dt::Real=1,
                      state0::AbstractVector{<:Real}=zeros(getDimension(domain)),
                      state1::AbstractVector{<:Real}=zeros(getDimension(domain)),
                      tspan::Tuple{<:Real,<:Real}=_calcTimespan(Xi),
                      problemType = DiscreteProblem )
  #
  datatuple = if 2 < length(Xi)
    datavec = getDimension.([_maketuplebeyond2args(Xi...)...]) .|> x->zeros(x)
    (data,datavec...,)
  else
    data
  end
  # forward time problem
  fproblem = problemType(f, state0, tspan, datatuple, dt=dt )
  # backward time problem
  bproblem = problemType(f, state1, (tspan[2],tspan[1]), datatuple, dt=-dt )
  # build the IIF recognizable object
  ODERelative( domain, fproblem, bproblem, datatuple, getSample )
end


ODERelative(dfg::AbstractDFG,
            labels::AbstractVector{Symbol},
            domain::Type{<:InferenceVariable},
            f::Function,
            data=()->();
            Xi::AbstractArray{<:DFGVariable}=getVariable.(dfg, labels),
            dt::Real=1,
            state0::AbstractVector{<:Real}=zeros(getDimension(domain)),
            state1::AbstractVector{<:Real}=zeros(getDimension(domain)),
            tspan::Tuple{<:Real,<:Real}=_calcTimespan(Xi),
            problemType = DiscreteProblem ) = ODERelative( Xi, domain, f, data; dt=dt, state0=state0, state1=state1, tspan=tspan, problemType=problemType  )
#
#


# Xtra splat are variable points (X3::Matrix, X4::Matrix,...)
function _solveFactorODE!(measArr, prob, u0pts, i, Xtra...)
  # should more variables be included in calculation
  for (xid, xtra) in enumerate(Xtra)
    # update the data register before ODE solver calls the function
    prob.p[xid+1][:] = Xtra[xid][:,i]
  end
  
  # set the initial condition
  prob.u0[:] = u0pts[:,i]
  sol = DifferentialEquations.solve(prob)
  
  # extract solution from solved ode
  measArr[:,i] = sol.u[end]
  sol
end

# FIXME see #1025, `multihypo=` will not work properly yet
function getSample( oder::ODERelative, 
                    N::Int=1, 
                    fmd_...)
  #

  # how many trajectories to propagate?
  meas = zeros(getDimension(fmd_[1].fullvariables[2]), N)
  
  # pick forward or backward direction
  prob = oder.forwardProblem
  # buffer manifold operations for use during factor evaluation
  addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks( getManifolds(fmd_[1].fullvariables[2]) )
  # set boundary condition
  u0pts = if fmd_[1].solvefor == DFG.getLabel(fmd_[1].fullvariables[1])
    # backward direction
    prob = oder.backwardProblem
    addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks( getManifolds(fmd_[1].fullvariables[1]) )
    getBelief( fmd_[1].fullvariables[2] ) |> getPoints
  else
    # forward backward
    getBelief( fmd_[1].fullvariables[1] ) |> getPoints
  end

  # solve likely elements
  for i in 1:N
    _solveFactorODE!(meas, prob, u0pts, i, _maketuplebeyond2args(fmd_[1].arrRef...)...)
  end

  return (meas, diffOp)
end
# getDimension(oderel.domain)


# FIXME see #1025, `multihypo=` will not work properly yet
function (oderel::ODERelative)( res::AbstractVector{<:Real},
                                fmd::FactorMetadata,
                                idx::Int,
                                meas::Tuple,
                                X...)
  #
  
  # work on-manifold
  diffOp = meas[2]
  # if backwardSolve else forward
  
  # check direction
  # TODO put solveforIdx in FMD?
  solveforIdx = 1
  if fmd.solvefor == DFG.getLabel(fmd.fullvariables[2])
    solveforIdx = 2
  elseif fmd.solvefor in _maketuplebeyond2args(fmd.variablelist...)
    # need to recalculate new ODE (forward) for change in parameters (solving for 3rd or higher variable)
    solveforIdx = 2
    # use forward solve for all solvefor not in [1;2]
    u0pts = getBelief(fmd.fullvariables[1]) |> getPoints
    # update parameters for additional variables
    _solveFactorODE!(meas[1], oderel.forwardProblem, u0pts, idx, _maketuplebeyond2args(X...)...)
  end

  # find the difference between measured and predicted.
  ## assuming the ODE integrated from current X1 through to predicted X2 (ie `meas[1][:,idx]`)
  ## FIXME, obviously this is not going to work for more compilcated groups/manifolds -- must fix this soon!
  for i in 1:size(X[2],1)
    # diffop( test, reference )   <===>   ΔX = test \ reference
    res[i] = diffOp[i]( X[solveforIdx][i,idx], meas[1][i,idx] )
  end
  res
end



## the function
# ode.problem.f.f


#