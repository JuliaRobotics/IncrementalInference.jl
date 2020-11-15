# Basic test of ODERelative


using Test
using DifferentialEquations
using IncrementalInference
using Dates
using Statistics

## plotting functions

using Plots
using Cairo, RoMEPlotting
Gadfly.set_default_plot_size(25cm,20cm)

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

doautoinit!(fg, :x0)

prev = :x0

for i in 1:3

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, ContinuousScalar, timestamp=DateTime(2000,1,1,0,0,5*i))
  oder = IIF.ODERelative( fg, [prev; nextSym], 
                          ContinuousEuclid{1}, 
                          firstOrder!,
                          tstForce,
                          dt=0.05, 
                          problemType=ODEProblem )
  #
  addFactor!( fg, [prev;nextSym], oder, graphinit=false )
  initManual!(fg, nextSym, zeros(1,100))

  prev = nextSym
end


## basic sample test

meas = freshSamples(fg, :x0x1f1, 10)
@test size(meas[1],1) == 1
@test size(meas[1],2) == 10


## do all forward solutions
pts, = freshSamples(fg, :x0f1, 100)
initManual!(fg, :x0, pts)
pts = approxConv(fg, :x0x1f1, :x1)
@test 0.3 < Statistics.mean(pts) < 0.4


# check that the reverse solve also works

initManual!(fg, :x1, pts)
pts = approxConv(fg, :x0x1f1, :x0)

# check the reverse solve to be relatively accurate
@test pts - (getBelief(fg, :x0) |> getPoints) |> norm < 1e-4


##

oder_ = ODERelative( fg, [:x0; :x3], 
                    ContinuousEuclid{1}, 
                    firstOrder!,
                    tstForce, 
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0]
sl = DifferentialEquations.solve(oder_.forwardProblem)

##


# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",layout=(1,1))

# for lb in [:x0; :x1;:x2;:x3]
#   x = getTimestamp(getVariable(fg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [0;1]
#   Plots.plot!(xx, yy, show=true)
# end


##


tfg = initfg()
pts = approxConv(fg, :x0f1, :x3, setPPE=true, tfg=tfg)
# initManual!(tfg, :x3, pts)


##


@test getPPE(tfg, :x0).suggested - sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(tfg, :x1).suggested - sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(tfg, :x2).suggested - sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test       Statistics.mean(pts) - sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix)[1] < 0.1


##

# plotKDE(tfg, [:x0;:x1;:x2;:x3])


## Now test a full solve

solveTree!(fg);


##


@test getPPE(fg, :x0).suggested - sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(fg, :x1).suggested - sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(fg, :x2).suggested - sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(fg, :x3).suggested - sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1


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
  ω = 0.7
  β = -0.3
  dstate[2] = β*state[2] - ω*state[1] + force(t)
  dstate[1] = state[2]
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0


## build a representative factor graph with ODE built inside

fg = initfg()

# the starting points and "0 seconds"
addVariable!(fg, :x0, ContinuousEuclid{2}, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(MvNormal([1;0],0.01*diagm(ones(2)))))



##

prev = :x0
DT = 2

for i in 1:7

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, ContinuousEuclid{2}, timestamp=DateTime(2000,1,1,0,0,DT*i))
  oder = ODERelative( fg, [prev; nextSym], 
                      ContinuousEuclid{2}, 
                      dampedOscillator!,
                      tstForce, 
                      # (state, var)->(state[1] = var[1]),
                      # (var, state)->(var[1] = state[1]),
                      dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev;nextSym], oder )

  prev = nextSym
end


##

tfg = initfg()
for s in ls(fg)
  initManual!(fg, s, zeros(2,100))
end

pts = approxConv(fg, :x0f1, :x7, setPPE=true, tfg=tfg)
# initManual!(tfg, :x7, pts)



##

plotKDE(tfg, ls(fg) |> sortDFG, dims=[1] )

##


oder_ = ODERelative( fg, [:x0; :x7], 
                    ContinuousEuclid{2}, 
                    dampedOscillator!,
                    tstForce, 
                    # (state, var)->(state[1] = var[1]),
                    # (var, state)->(var[1] = state[1]),
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0;0.0]
sl = DifferentialEquations.solve(oder_.forwardProblem)



##



Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))

for lb in sortDFG(ls(fg))
  x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
  xx = [x;x]
  yy = [-1;1]
  Plots.plot!(xx, yy, show=true)
end


##

end





##

@testset "Parameterized Damped Oscillator ODERelative" begin

## setup some example dynamics

# Lets build an damped oscillator to demonstrate the process in state space
# https://en.wikipedia.org/wiki/Harmonic_oscillator
# ddx/ddt = β dx/dt  -  ω x  +  force[t]
# dx/dt   = dx/dt
# force_ωβ = (data, ωβ)
function dampedOscillatorParametrized!(dstate, state, force_ωβ, t)
  # 3rd variable in this factor graph test example
  force = force_ωβ[1]
  @show ω     = force_ωβ[2][1]
  @show β     = force_ωβ[2][2]
  # classic ODE between first and second fg variables
  dstate[2] = β*state[2] - ω*state[1] + force(t)
  dstate[1] = state[2]
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0


## build a representative factor graph with ODE built inside

fg = initfg()

# the starting points and "0 seconds"
addVariable!(fg, :x0, ContinuousEuclid{2}, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(MvNormal([1;0],0.01*diagm(ones(2)))))
doautoinit!(fg, :x0)

# and the new parameterized variable
ω = 0.7
β = -0.3

# these are the stochastic parameters
addVariable!(fg, :ωβ, ContinuousEuclid{2}) # timestamp should not matter
# pin with a simple prior
addFactor!(fg, [:ωβ], Prior(MvNormal([ω;β],0.0001*diagm(ones(2)))))
doautoinit!(fg, :ωβ)


##

prev = :x0
DT = 2

for i in 1:7

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, ContinuousEuclid{2}, timestamp=DateTime(2000,1,1,0,0,DT*i))
  oder = ODERelative( fg, [prev; nextSym; :ωβ], 
                      ContinuousEuclid{2}, 
                      dampedOscillatorParametrized!,
                      tstForce, # this is passed in as `force_ωβ[1]`
                      # (state, var)->(state[1] = var[1]),
                      # (var, state)->(var[1] = state[1]),
                      # dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev; nextSym; :ωβ], oder, graphinit=false )

  prev = nextSym
end


## test basic single approxConv

pts = approxConv(fg, :x0f1, :x1)


##

tfg = initfg()
for s in ls(fg)
  initManual!(fg, s, zeros(2,100))
end

# must initialize the parameters
pts = approxConv(fg, :ωβf1, :ωβ)
initManual!(fg, :ωβ, pts)

# project forward
forcepath = [:x0f1;]
push!(forcepath, :x0) 
push!(forcepath, :x0x1ωβf1) 
push!(forcepath, :x1)
push!(forcepath, :x1x2ωβf1)
push!(forcepath, :x2)
push!(forcepath, :x2x3ωβf1)
push!(forcepath, :x3)
push!(forcepath, :x3x4ωβf1)
push!(forcepath, :x4)
push!(forcepath, :x4x5ωβf1)
push!(forcepath, :x5)
push!(forcepath, :x5x6ωβf1)
push!(forcepath, :x6)
push!(forcepath, :x6x7ωβf1)
push!(forcepath, :x7)
pts = approxConv(fg, :x0f1, :x7, setPPE=true, tfg=tfg, path=forcepath)


##

plotKDE(tfg, ls(tfg) |> sortDFG, dims=[1] )


##

# getBelief(fg, :ωβ) |> getPoints

plotKDE(tfg, :ωβ)

##


oder_ = ODERelative( fg, [:x0; :x7; :ωβ], 
                    ContinuousEuclid{2}, 
                    dampedOscillatorParametrized!,
                    tstForce,
                    # (state, var)->(state[1] = var[1]),
                    # (var, state)->(var[1] = state[1]),
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0;0.0]
oder_.data[2] .= [ω;β]
sl = DifferentialEquations.solve(oder_.forwardProblem)



##



Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))

for lb in sortDFG(ls(fg))
  x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
  xx = [x;x]
  yy = [-1;1]
  Plots.plot!(xx, yy, show=true)
end


##

end





@error "ODERelative not tested for `multihypo=` case yet, see issue #1025"




#