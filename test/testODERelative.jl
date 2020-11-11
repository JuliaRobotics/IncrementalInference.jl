# Basic test of ODERelative


using Test
using DifferentialEquations
using IncrementalInference
using Dates
using Statistics

##

@testset "First order ODERelative" begin

##

function firstOrder!(dstate, state, force, t)
  β = -0.2
  dstate[1] = β*state[1] + force(t)
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0

## build a representative factor graph with ODE built inside

fg = initfg()
# the starting points and "0 seconds"
addVariable!(fg, :x0, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(Normal(1,0.01)))


prev = :x0

for i in 1:3

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,5*i))
  oder = ODERelative( fg, [prev; nextSym], 
                      ContinuousEuclid{1}, 
                      firstOrder!,
                      tstForce,
                      dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev;nextSym], oder )

  prev = nextSym
end


##

oder_ = ODERelative( fg, [:x0; :x3], 
                    ContinuousEuclid{1}, 
                    firstOrder!,
                    tstForce, 
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.problem.u0 .= [1.0]
sl = DifferentialEquations.solve(oder_.problem)

##

# using Plots
# using Cairo, RoMEPlotting
# Gadfly.set_default_plot_size(35cm,20cm)

##

# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]"],layout=(1,1))

# for lb in [:x1;:x2;:x3]
#   x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [0;1]
#   Plots.plot!(xx, yy, show=true)
# end

##


tfg = initfg()
pts = approxConv(fg, :x0f1, :x3, tfg=tfg)
initManual!(tfg, :x3, pts)
calcPPE(getVariable(tfg, :x3))



@test getPPE(tfg, :x0).suggested - sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(tfg, :x1).suggested - sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(tfg, :x2).suggested - sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test Statistics.mean(pts) - sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix)[1] < 0.1


##

# plotKDE(tfg, [:x0;:x1;:x2;:x3])

##

end




##

@testset "Damped Oscillator ODERelative" begin

## setup some example dynamics

# Lets build an damped oscillator to demonstrate the process in state space
# https://en.wikipedia.org/wiki/Harmonic_oscillator
# ddx/ddt = β dx/dt  -  ω x  +  force[t]
# dx/dt   = dx/dt
function dampedOscillator!(dstate, state, force, t)
  ω = 0.5
  β = -0.3
  dstate[1] = β*state[1] - ω*state[2] + force(t)
  dstate[2] = state[1]
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0


## build a representative factor graph with ODE built inside

fg = initfg()

# the starting points and "0 seconds"
addVariable!(fg, :x0, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(Normal(1,0.01)))



##

prev = :x0

for i in 1:3

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,5*i))
  oder = ODERelative( fg, [prev; nextSym], 
                      ContinuousEuclid{2}, 
                      dampedOscillator!,
                      tstForce, 
                      (state, var)->(state[2] = var[1]),
                      (var, state)->(var[1] = state[2]),
                      dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev;nextSym], oder )

  prev = nextSym
end


##

tfg = initfg()
pts = approxConv(fg, :x0f1, :x3, tfg=tfg)
initManual!(tfg, :x3, pts)


##

using Plots
using Cairo, RoMEPlotting
Gadfly.set_default_plot_size(35cm,20cm)

##

plotKDE(tfg, [:x0;:x1;:x2;:x3])

##



##

oder_ = ODERelative( fg, [:x0; :x3], 
                    ContinuousEuclid{2}, 
                    dampedOscillator!,
                    tstForce, 
                    (state, var)->(state[2] = var[1]),
                    (var, state)->(var[1] = state[2]),
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.problem.u0 .= [0.0;1.0]
sl = DifferentialEquations.solve(oder_.problem)

##


##

# using TesnsorCast

# sl.u

Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))

for lb in [:x1;:x2;:x3]
  x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
  xx = [x;x]
  yy = [-1;1]
  Plots.plot!(xx, yy, show=true)
end


##

end