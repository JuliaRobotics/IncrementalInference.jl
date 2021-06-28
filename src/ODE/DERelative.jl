
@info "IncrementalInference.jl is loading tools related to DifferentialEquations.jl"

using .DifferentialEquations

import .DifferentialEquations: solve

import IncrementalInference: getSample

export DERelative


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
struct DERelative{T <:InferenceVariable, P, D} <: AbstractRelativeRoots
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

# performance helper function, FIXME not compatible with all multihypo cases
_maketuplebeyond2args = (w1=nothing,w2=nothing,w3_...) -> (w3_...,)


function DERelative( Xi::AbstractVector{<:DFGVariable},
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
  DERelative( domain, fproblem, bproblem, datatuple, getSample )
end


DERelative( dfg::AbstractDFG,
            labels::AbstractVector{Symbol},
            domain::Type{<:InferenceVariable},
            f::Function,
            data=()->();
            Xi::AbstractArray{<:DFGVariable}=getVariable.(dfg, labels),
            dt::Real=1,
            state0::AbstractVector{<:Real}=zeros(getDimension(domain)),
            state1::AbstractVector{<:Real}=zeros(getDimension(domain)),
            tspan::Tuple{<:Real,<:Real}=_calcTimespan(Xi),
            problemType = DiscreteProblem ) = DERelative( Xi, domain, f, data; dt=dt, state0=state0, state1=state1, tspan=tspan, problemType=problemType  )
#
#


# Xtra splat are variable points (X3::Matrix, X4::Matrix,...)
function _solveFactorODE!(measArr, prob, u0pts, Xtra...)
  # should more variables be included in calculation
  for (xid, xtra) in enumerate(Xtra)
    # update the data register before ODE solver calls the function
    prob.p[xid+1][:] = Xtra[xid][:]
  end
  
  # set the initial condition
  prob.u0[:] = u0pts[:]
  sol = DifferentialEquations.solve(prob)
  
  # extract solution from solved ode
  measArr[:] = sol.u[end]
  sol
end

# FIXME see #1025, `multihypo=` will not work properly yet
function getSample( cf::CalcFactor{<:DERelative}, 
                    N::Int=1 )
  #

  oder = cf.factor
  fmd_ = cf.metadata

  # how many trajectories to propagate?
  # @show getLabel(fmd_.fullvariables[2]), getDimension(fmd_.fullvariables[2])
  meas = [zeros(getDimension(fmd_.fullvariables[2])) for _ in 1:N]
  
  # pick forward or backward direction
  prob = oder.forwardProblem
  # buffer manifold operations for use during factor evaluation
  addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks( convert(Tuple, getManifold(getVariableType(fmd_.fullvariables[2]))) )
  
  # set boundary condition
  u0pts = if fmd_.solvefor == DFG.getLabel(fmd_.fullvariables[1])
    # backward direction
    prob = oder.backwardProblem
    addOp, diffOp, _, _ = AMP.buildHybridManifoldCallbacks( convert(Tuple, getManifold(getVariableType(fmd_.fullvariables[1]))) )
    getBelief( fmd_.fullvariables[2] ) |> getPoints
  else
    # forward backward
    getBelief( fmd_.fullvariables[1] ) |> getPoints
  end

  # solve likely elements
  for i in 1:N
    idxArr = getindex.(fmd_.arrRef, i)
    _solveFactorODE!(meas[i], prob, u0pts[i], _maketuplebeyond2args(idxArr...)...)
    # _solveFactorODE!(meas, prob, u0pts, i, _maketuplebeyond2args(fmd_.arrRef...)...)
  end

  @info "this is diffOp" diffOp
  return (meas, diffOp)
end
# getDimension(oderel.domain)


# NOTE see #1025, CalcFactor should fix `multihypo=` in `cf.metadata` fields
function (cf::CalcFactor{<:DERelative})(meas1,
                                        diffOp,
                                        X...)
  #
  
  oderel = cf.factor
  
  # work on-manifold
  # diffOp = meas[2]
  # if backwardSolve else forward
  
  # check direction
  # TODO put solveforIdx in FMD?
  solveforIdx = 1
  if cf.metadata.solvefor == DFG.getLabel(cf.metadata.fullvariables[2])
    solveforIdx = 2
  elseif cf.metadata.solvefor in _maketuplebeyond2args(cf.metadata.variablelist...)
    # need to recalculate new ODE (forward) for change in parameters (solving for 3rd or higher variable)
    solveforIdx = 2
    # use forward solve for all solvefor not in [1;2]
    u0pts = getBelief(cf.metadata.fullvariables[1]) |> getPoints
    # update parameters for additional variables
    _solveFactorODE!(meas1, oderel.forwardProblem, u0pts[cf._sampleIdx], _maketuplebeyond2args(X...)...)
  end

  # find the difference between measured and predicted.
  ## assuming the ODE integrated from current X1 through to predicted X2 (ie `meas1[:,idx]`)
  ## FIXME, obviously this is not going to work for more compilcated groups/manifolds -- must fix this soon!
  # @show cf._sampleIdx, solveforIdx, meas1

  #FIXME 
  res = zeros(size(X[2],1))
  for i in 1:size(X[2],1)
    # diffop( test, reference )   <===>   ΔX = test \ reference
    @info "what what ode" i solveforIdx diffOp X[solveforIdx][i] meas1[i]
    res[i] = diffOp( X[solveforIdx][i], meas1[i] )
  end
  res
end



## the function
# ode.problem.f.f


#