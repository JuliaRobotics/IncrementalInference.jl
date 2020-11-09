
@info "Loading DifferentialEquations.jl tools in IncrementalInference.jl"

using .DifferentialEquations

import IncrementalInference: getSample

export ODERelative


struct ODERelative{T <:InferenceVariable, P} <: AbstractRelativeRoots
  domain::Type{T}
  problem::P
  specialSampler::Function
end



"""
    $SIGNATURES

Calculate a DifferentialEquations.jl ready `tspan::Tuple{Float64,Float64}` from DFGVariables.

Notes
- Can change numerical data return type using an additional first argument, `_calcTimespan(Float32, Xi)`.

DevNotes
- # TODO does not yet incorporate Xi.nanosecond field.
"""
function _calcTimespan( elemType::Type{<:Real},
                        Xi::AbstractVector{<:DFGVariable} )
  #
  tsmps = getTimestamp.(Xi)
  toffs = (tsmps .- tsmps[1]) .|> x-> elemType(x.value*1e-3) 
  return (toffs...,)
end

_calcTimespan(Xi::AbstractVector{<:DFGVariable}) = _calcTimespan(Float64, Xi)


function ODERelative( Xi::AbstractVector{<:DFGVariable},
                      domain::Type{<:InferenceVariable},
                      f::Function;
                      dt::Real=1,
                      u0::AbstractVector{<:Real}=zeros(getDimension(domain)),
                      tspan::Tuple{<:Real,<:Real}=_calcTimespan(Xi)  )
  #
  prob = DiscreteProblem(f, u0, tspan, dt=dt)
  ODERelative( domain, prob, getSample )
end

ODERelative(dfg::AbstractDFG,
            labels::AbstractVector{Symbol},
            domain::Type{<:InferenceVariable},
            f::Function;
            Xi::AbstractArray{<:DFGVariable}=getVariable.(dfg, labels),
            dt::Real=1,
            u0::AbstractVector{<:Real}=zeros(getDimension(domain)),
            tspan::Tuple{<:Real,<:Real}=_calcTimespan(Xi)  ) = ODERelative( Xi, domain, f; dt=dt, u0=u0, tspan=tspan  )
#

function getSample( oder::ODERelative, 
                    N::Int=1, 
                    fmd_...)
  #
  # how many trajectories to propagate?
  meas = zeros(getDimension(oder.domain),N)
  
  u0pts = getBelief(fmd_[1].fullvariables[1]) |> getPoints

  # solve likely elements
  for i in 1:N
    # set the initial condition
    oder.problem.u0 .= u0pts[:,i]

    # solve the ode
    meas[:,i] = DifferentialEquations.solve(oder.problem).u[end]
  end

  return (meas,)
end


function (oderel::ODERelative)( res::AbstractVector{<:Real},
                                fmd::FactorMetadata,
                                idx::Int,
                                meas::Tuple,
                                X...)
  #
  fill!(res, 0.0)
  res .= meas[1][:,idx] - X[2][:,idx]
end

## the function
# ode.problem.f.f


#