



function SDERelative(
  Xi::AbstractVector{<:DFGVariable},
  domain::Type{<:InferenceVariable},
  f::Function, # imuKinematicQ!
  fσ::Function, # sigma_imuKinematicQ!
  data = () -> ();
  dt::Real = 1,
  state0::AbstractVector{<:Real} = zeros(getDimension(domain)),
  state1::AbstractVector{<:Real} = zeros(getDimension(domain)),
  tspan::Tuple{<:Real, <:Real} = _calcTimespan(Xi),
  problemType = DiscreteProblem,
)
  #
  datatuple = if 2 < length(Xi)
    datavec = getDimension.([_maketuplebeyond2args(Xi...)...]) .|> x -> zeros(x)
    (data, datavec...)
  else
    data
  end

  prob = SDEProblem(f, fσ, u0, tspan, Ref((gyro=gyros_t, accel=accels_t)))

  # forward time problem
  fproblem = problemType(f, state0, tspan, datatuple; dt = dt)
  # backward time problem
  bproblem = problemType(f, state1, (tspan[2], tspan[1]), datatuple; dt = -dt)
  # build the IIF recognizable object
  return DERelative(domain, fproblem, bproblem, datatuple, getSample)
end