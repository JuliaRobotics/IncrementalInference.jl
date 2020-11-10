# Basic test of ODERelative


using Test
using DifferentialEquations
using IncrementalInference
using Dates


@testset "Basic test of ODERelative" begin

fg = initfg()

# the starting points and "0 seconds"
addVariable!(fg, :x0, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(Normal(1,0.1)))

# another point in the trajectory 10 seconds later
addVariable!(fg, :x1, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,5))

# another point in the trajectory 10 seconds later
addVariable!(fg, :x2, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,10))

# Lets build an damped oscillator to demonstrate the process in state space
# https://en.wikipedia.org/wiki/Harmonic_oscillator
# ddx/ddt = β dx/dt  -  ω x  +  force[t]
# dx/dt   = dx/dt
function dampedOscillator!(dstate, state, force, t)
  ω = 3
  β = -0.5
  dstate[1] = β*state[1] - ω*state[2] + force(t)
  dstate[2] = state[1]
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0


oder1 = ODERelative( fg, [:x0; :x1], 
                    ContinuousEuclid{2}, 
                    dampedOscillator!,
                    tstForce, 
                    (state, var)->(state[2] = var[1]),
                    (var, state)->(var[1] = state[2]),
                    dt=0.05, 
                    problemType=ODEProblem )
#

oder2 = ODERelative( fg, [:x1; :x2], 
                    ContinuousEuclid{2}, 
                    dampedOscillator!,
                    tstForce, 
                    (state, var)->(state[2] = var[1]),
                    (var, state)->(var[1] = state[2]),
                    dt=0.05, 
                    problemType=ODEProblem )


addFactor!( fg, [:x0;:x1], oder1 )

freshSamples(fg, :x0x1f1, 1)


# doautoinit!(fg, :x1)

addFactor!( fg, [:x1; :x2], oder2 )

# doautoinit!(fg, :x2)



oder1.problem.u0 .= [0.0;1.0]
sl = DifferentialEquations.solve(oder1.problem)


using Plots
Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))


end