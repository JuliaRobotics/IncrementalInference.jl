# Basic test of DERelative


using Test
using DifferentialEquations
using IncrementalInference
using Dates
using Statistics
using TensorCast

## plotting functions

# using Plots
# using Cairo, RoMEPlotting
# Gadfly.set_default_plot_size(25cm,20cm)

##

@testset "First order DERelative" begin

##

# a user specified ODE in standard form
# inplace `xdot = f(x, u, t)`
# if linear, `xdot = F*x(t) + G*u(t)`
function firstOrder!(dstate, state, u, t)
  β = -0.2
  dstate[1] = β*state[1] + u(t)
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0


## build a representative factor graph with ODE built inside

fg = initfg()
# the starting points and "0 seconds"
# `accurate_time = trunc(getDatetime(var), Second) + (1e-9*getNstime(var) % 1)`
addVariable!(fg, :x0, Position{1}, timestamp=DateTime(2000,1,1,0,0,0)) 
# pin with a simple prior
addFactor!(fg, [:x0], Prior(Normal(1,0.01)))

doautoinit!(fg, :x0)

prev = :x0

for i in 1:3

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, Position{1}, timestamp=DateTime(2000,1,1,0,0,5*i))
  # build factor against manifold Manifolds.TranslationGroup(1)
  ode_fac = IIF.DERelative(
    fg, [prev; nextSym],
    Position{1},
    firstOrder!,
    tstForce,
    dt=0.05,
    problemType=ODEProblem,
    keepSolution = i == 1,
  )
  #
  addFactor!( fg, [prev;nextSym], ode_fac, graphinit=false, keepCalcFactor = i == 1 )
  initVariable!(fg, nextSym, [0.1*randn(1) for _ in 1:100])

  prev = nextSym
end



## raw test against DiffEq API directly

oder_ = DERelative( fg, [:x0; :x3], 
                    Position{1}, 
                    firstOrder!,
                    tstForce, 
                    dt=0.05, 
                    problemType=ODEProblem )



##

oder_.forwardProblem.u0 .= [1.0]
sl = DifferentialEquations.solve(oder_.forwardProblem)

x0_val_ref = sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix)
x1_val_ref = sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix)
x2_val_ref = sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix)
x3_val_ref = sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix)


# sanity check that the backward problem runs through
oder_.backwardProblem.u0 .= [0.049788]
retcode = DifferentialEquations.solve(oder_.backwardProblem)

@test isapprox(1, retcode.u[end][1]; atol=1e-4)


## one layer wrapped API test through IIFExt to DiffEq

f = getFactor(fg, intersect(ls(fg,:x0),ls(fg,:x1))[1] )
fc = getObservation(f)

fprob = fc.forwardProblem

meas = zeros(getDimension(getVariable(fg, :x1)))
u0pts = getPoints(getBelief(fg, :x0))[1]
res = IncrementalInference._solveFactorODE!(meas, fprob, u0pts)

@test isapprox( 5, res.t[end]-res.t[1]; atol=1e-6)
@test isapprox( x0_val_ref, res.u[1]; atol=0.1)
@test isapprox( x1_val_ref, res.u[end]; atol=0.1)


## basic sample test

meas = sampleFactor(fg, :x0x1f1, 10)
@test size(meas[1][1],1) == 1
@test size(meas,1) == 10


## do all forward solutions
chnl = Channel(100)


pts = sampleFactor(fg, :x0f1, 100)

initVariable!(fg, :x0, pts)
pts_ = approxConv(fg, :x0x1f1, :x1; keepCalcFactor=chnl)
@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.mean(pts) < 0.4

ccf = take!(chnl)
@test 50 < length(ccf.factor.keepSolution)

## check that the reverse solve also works

initVariable!(fg, :x1, pts_)
pts_ = approxConv(fg, :x0x1f1, :x0)
@cast pts[i,j] := pts_[j][i]

# check the reverse solve to be relatively accurate
ref_ = (getBelief(fg, :x0) |> getPoints)
@cast ref[i,j] := ref_[j][i]
@test norm(pts - ref) < 1e-4


##

# use Makie instead
# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",layout=(1,1))

# for lb in [:x0; :x1;:x2;:x3]
#   x = getTimestamp(getVariable(fg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [0;1]
#   Plots.plot!(xx, yy, show=true)
# end


## temp graph solve check

tfg  = initfg()
tx3_ = approxConvBelief(fg, :x0f1, :x3; tfg)
pts_ = getPoints(tx3_)
# initVariable!(tfg, :x3, pts)

@cast pts[i,j] := pts_[j][i]

@test isapprox( x0_val_ref, calcMeanMaxSuggested(tfg, :x0, :default).suggested ; atol = 0.1)
@test isapprox( x1_val_ref, calcMeanMaxSuggested(tfg, :x1, :default).suggested ; atol = 0.1)
@test isapprox( x2_val_ref, calcMeanMaxSuggested(tfg, :x2, :default).suggested ; atol = 0.1)
@test isapprox( x3_val_ref, mean(tx3_); atol=0.1)

# using KernelDensityEstimatePlotting
# plotKDE(tfg, [:x0;:x1;:x2;:x3])


## check if variables are initialized (only works for graphinit)

@test isInitialized(fg, :x0)
@test isapprox( x0_val_ref, mean(getBelief(fg[:x0])); atol=0.1)

@test isInitialized(fg, :x1)
# @test isapprox( x1_val_ref, mean(getBelief(fg[:x1])); atol=0.1)

X2_ = approxConvBelief(fg, :x1x2f1, :x2)
@test isapprox( x2_val_ref, mean(X2_); atol=0.1)

# FIXME, X2 and X3 are wrongly initialized to zero above
# X2_ = approxConvBelief(fg, :x2x3f1, :x2)
# @test isapprox( x2_val_ref, mean(X2_); atol=0.1)
# @enter approxConvBelief(fg, :x2x3f1, :x2)

factors = getFactor.(fg, IIF.listNeighbors(fg, :x2))
dens = ManifoldKernelDensity[]
ipc = IIF.proposalbeliefs!(fg, :x2, factors, dens)

# 
mkd = *(dens...)

@test isapprox( x2_val_ref, mean(mkd); atol=0.1)

X2_,_ = propagateBelief(fg, :x2, :)
@test isapprox( x2_val_ref, mean(X2_); atol=0.1)
# @enter propagateBelief(fg, :x2, :)

@test isInitialized(fg, :x2)
@test isInitialized(fg, :x3)

# FIXME, wrongly initialized X2 and X3 to near zero above
# @test isapprox( x2_val_ref, mean(getBelief(fg[:x2])); atol=0.1)
# @test isapprox( x3_val_ref, mean(getBelief(fg[:x3])); atol=0.1) # happens to be near zero


## Now test a full graph solve

smtasks = Task[]
tree = solveTree!(fg; smtasks, recordcliqs=ls(fg));

oder_.keepSolution


hists = fetchCliqHistoryAll!(smtasks)

printCSMHistoryLogical(hists)


##

# intended steps at writing are 5, 6 (upsolve) 
_, csmc = repeatCSMStep!(hists[1], 5; duplicate=true)
@test isapprox( 1, calcMeanMaxSuggested(csmc.cliqSubFg, :x0, :default).suggested[1]; atol=0.1 )
nval_x0 = mean(getBelief(csmc.cliqSubFg, :x0))
@test isapprox( x0_val_ref, nval_x0; atol=0.1 )

nval_x1 = mean(getBelief(csmc.cliqSubFg, :x1))
@test isapprox( x1_val_ref, nval_x1; atol=0.1 )


sfg = deepcopy( hists[1][6][4].cliqSubFg )
dens, ipc = propagateBelief( sfg,  :x0,  :;)
@test isapprox( x0_val_ref, mean(dens); atol=0.1)

@test isapprox( x0_val_ref, mean(getBelief(sfg[:x0])); atol=0.1)
# @test isapprox( x2_val_ref, mean(getBelief(sfg[:x2])); atol=0.1) # TODO DELETE THIS LINE

dens, ipc = propagateBelief( sfg,  :x1,  :;)
@test isapprox( x1_val_ref, mean(dens); atol=0.1)
# @enter propagateBelief(sfg,  :x1,  :)

_, csmc = repeatCSMStep!(hists[1], 6; duplicate=true)
# @enter repeatCSMStep!(hists[1], 6; duplicate=true)
@test isapprox( x0_val_ref, calcMeanMaxSuggested(csmc.cliqSubFg, :x0, :default).suggested; atol=0.1 )
nval_x0 = mean(getBelief(csmc.cliqSubFg, :x0))
@test isapprox( x0_val_ref, nval_x0; atol=0.1 )

nval_x0 = mean(getBelief(csmc.cliqSubFg, :x0))
@test isapprox( x0_val_ref, nval_x0; atol=0.1 )


# TODO CHECK vnd.val points istype SArray???

# intended steps at writing are 11,12 (post-root clique downsolve)
val0 = calcMeanMaxSuggested(hists[1][11][4].cliqSubFg, :x0, :default).suggested
@test isapprox( x0_val_ref, val0; atol=0.1)
val0 = calcMeanMaxSuggested(hists[1][12][4].cliqSubFg, :x0, :default).suggested
@test isapprox( x0_val_ref, val0; atol=0.1)


##

@test isapprox( calcMeanMaxSuggested(fg, :x0, :default).suggested, x0_val_ref; atol = 0.1)
@test isapprox( calcMeanMaxSuggested(fg, :x1, :default).suggested, x1_val_ref; atol = 0.1)
@test isapprox( calcMeanMaxSuggested(fg, :x2, :default).suggested, x2_val_ref; atol = 0.1)
@test isapprox( calcMeanMaxSuggested(fg, :x3, :default).suggested, x3_val_ref; atol = 0.1)

##

end


##

@testset "Damped Oscillator DERelative" begin

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
addVariable!(fg, :x0, Position{2}, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(MvNormal([1;0],0.01*diagm(ones(2)))))


##

prev = :x0
DT = 2

for i in 1:7

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, Position{2}, timestamp=DateTime(2000,1,1,0,0,DT*i))
  oder = DERelative( fg, [prev; nextSym], 
                      Position{2}, 
                      dampedOscillator!,
                      tstForce, 
                      # (state, var)->(state[1] = var[1]),
                      # (var, state)->(var[1] = state[1]),
                      dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev;nextSym], oder; graphinit=false )

  prev = nextSym
end


##

oder_ = DERelative( fg, [:x0; :x7], 
                    Position{2}, 
                    dampedOscillator!,
                    tstForce, 
                    # (state, var)->(state[1] = var[1]),
                    # (var, state)->(var[1] = state[1]),
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0;0.0]
sl = DifferentialEquations.solve(oder_.forwardProblem)

## Initialize the rest of the variables

initAll!(fg)

## check the solve values are correct

x0_val_ref = sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix)
x1_val_ref = sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix)
x2_val_ref = sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix)
x3_val_ref = sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix)
x4_val_ref = sl(getVariable(fg, :x4) |> getTimestamp |> DateTime |> datetime2unix)
x5_val_ref = sl(getVariable(fg, :x5) |> getTimestamp |> DateTime |> datetime2unix)
x6_val_ref = sl(getVariable(fg, :x6) |> getTimestamp |> DateTime |> datetime2unix)
x7_val_ref = sl(getVariable(fg, :x7) |> getTimestamp |> DateTime |> datetime2unix)

##

@test isapprox( calcMeanMaxSuggested(fg, :x0, :default).suggested, x0_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x1, :default).suggested, x1_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x2, :default).suggested, x2_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x3, :default).suggested, x3_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x4, :default).suggested, x4_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x5, :default).suggested, x5_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x6, :default).suggested, x6_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x7, :default).suggested, x7_val_ref; atol=0.2)


## check forward and backward solving

pts_ = approxConv(fg, :x0f1, :x0)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [1;0]) < 0.3

initVariable!(fg, :x0, pts_)
X0_ = deepcopy(pts)

pts_ = approxConv(fg, :x0x1f1, :x1)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [0;-0.6]) < 0.4

# now check the reverse direction solving
initVariable!(fg, :x1, pts_)
pts_ = approxConv(fg, :x0x1f1, :x0)
@cast pts[i,j] := pts_[j][i]

# check forward then backward convolves are reversible
@test isapprox(0, norm(X0_ - pts); atol=1e-2)


##

# plotKDE(tfg, ls(fg) |> sortDFG, dims=[1] )

##

tfg = initfg()
# for s in ls(fg)
#   initVariable!(fg, s, [0.1.*zeros(2) for _ in 1:100])
# end

pts = approxConv(fg, :x0f1, :x7, tfg=tfg)
initVariable!(tfg, :x7, pts)

##

@test isapprox( calcMeanMaxSuggested(tfg, :x0, :default).suggested, x0_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(tfg, :x1, :default).suggested, x1_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(tfg, :x2, :default).suggested, x2_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(tfg, :x3, :default).suggested, x3_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(tfg, :x4, :default).suggested, x4_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(tfg, :x5, :default).suggested, x5_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(tfg, :x6, :default).suggested, x6_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(tfg, :x7, :default).suggested, x7_val_ref; atol=0.2)

##

@error "Disabling useMsgLikelihood for DERelative test, follow fix on #1010 as rough guide"
getSolverParams(fg).useMsgLikelihoods = false

smtasks = Task[]
tree = solveTree!(fg; recordcliqs=ls(fg), smtasks);

hists = fetchCliqHistoryAll!(smtasks)
printCSMHistoryLogical(hists)

_, csmc = repeatCSMStep!(hists[2], 6; duplicate=true);


## 

# solveTree has weird problem in breaking correct init and inserting zeros???
@test isapprox( calcMeanMaxSuggested(fg, :x0, :default).suggested, x0_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x1, :default).suggested, x1_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x2, :default).suggested, x2_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x3, :default).suggested, x3_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x4, :default).suggested, x4_val_ref; atol=0.2)

@test isapprox( calcMeanMaxSuggested(fg, :x5, :default).suggested, x5_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x6, :default).suggested, x6_val_ref; atol=0.2)
@test isapprox( calcMeanMaxSuggested(fg, :x7, :default).suggested, x7_val_ref; atol=0.2)


##

# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))

# for lb in sortDFG(ls(fg))
#   x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [-1;1]
#   Plots.plot!(xx, yy, show=true)
# end

##

end





##

@testset "Parameterized Damped Oscillator DERelative (n-ary factor)" begin

## setup some example dynamics

# Lets build an damped oscillator to demonstrate the process in state space
# https://en.wikipedia.org/wiki/Harmonic_oscillator
# ddx/ddt = β dx/dt  -  ω x  +  force[t]
# dx/dt   = dx/dt
# force_ωβ = (data, ωβ)
function dampedOscillatorParametrized!(dstate, state, force_ωβ, t)
  # 3rd variable in this factor graph test example
  force = force_ωβ[1]
  ω     = force_ωβ[2][1]
  β     = force_ωβ[2][2]
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
addVariable!(fg, :x0, Position{2}, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(MvNormal([1;0],0.01*diagm(ones(2)))))
doautoinit!(fg, :x0)

# and the new parameterized variable
ω = 0.7
β = -0.3

# these are the stochastic parameters
addVariable!(fg, :ωβ, Position{2}) # timestamp should not matter
# pin with a simple prior
addFactor!(fg, [:ωβ], Prior(MvNormal([ω;β],0.0001*diagm(ones(2)))))
doautoinit!(fg, :ωβ)


##

prev = :x0
DT = 2

for i in 1:7

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, Position{2}, timestamp=DateTime(2000,1,1,0,0,DT*i))
  oder = DERelative( fg, [prev; nextSym; :ωβ], 
                      Position{2}, 
                      dampedOscillatorParametrized!,
                      tstForce, # this is passed in as `force_ωβ[1]`
                      # (state, var)->(state[1] = var[1]),
                      # (var, state)->(var[1] = state[1]),
                      # dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev; nextSym; :ωβ], oder, graphinit=false, inflation=0.01 )

  prev = nextSym
end


## check forward and backward solving

pts_ = approxConv(fg, :x0f1, :x0)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [1;0]) < 0.3

initVariable!(fg, :x0, pts_)
X0_ = deepcopy(pts)

pts_ = approxConv(fg, :x0x1ωβf1, :x1)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [0;-0.6]) < 0.4

# now check the reverse direction solving
initVariable!(fg, :x1, pts_)

# failing here
pts_ = approxConv(fg, :x0x1ωβf1, :x0)
@cast pts[i,j] := pts_[j][i]

@test (X0_ - pts) |> norm < 1e-2


##

tfg = initfg()
# for s in ls(fg)
#   initVariable!(fg, s, [zeros(2) for _ in 1:100])
# end

# must initialize the parameters
pts = approxConv(fg, :ωβf1, :ωβ)
initVariable!(fg, :ωβ, pts)

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
pts = approxConv(fg, :x0f1, :x7, tfg=tfg, path=forcepath)


##

# plotKDE(tfg, ls(tfg) |> sortDFG, dims=[1] )


##

# getBelief(fg, :ωβ) |> getPoints

# plotKDE(tfg, :ωβ)

##


oder_ = DERelative( fg, [:x0; :x7; :ωβ], 
                    Position{2}, 
                    dampedOscillatorParametrized!,
                    tstForce,
                    # (state, var)->(state[1] = var[1]),
                    # (var, state)->(var[1] = state[1]),
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0;0.0]
oder_.data[2] .= [ω;β]
sl = DifferentialEquations.solve(oder_.forwardProblem)



## check the approxConv is working right

for sym in setdiff(ls(tfg), [:ωβ])
  @test calcMeanMaxSuggested(tfg, sym, :default).suggested - sl(getVariable(fg, sym) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.2
end


## 


# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))

# for lb in sortDFG(ls(fg))
#   x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [-1;1]
#   Plots.plot!(xx, yy, show=true)
# end


## test convolution to the parameter (third) variable

# easy test with good starting points
pts = approxConv(fg, :ωβf1, :ωβ)
initVariable!(fg, :ωβ, pts)

# make sure the other variables are in the right place
pts_ = getBelief(fg, :x0) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test Statistics.mean(pts, dims=2) - [1;0] |> norm < 0.1

pts_ = getBelief(fg, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test Statistics.mean(pts, dims=2) - [0;-0.6] |> norm < 0.2


pts_ = approxConv(fg, :x0x1ωβf1, :ωβ)
@cast pts[i,j] := pts_[j][i]
@test Statistics.mean(pts, dims=2) - [0.7;-0.3] |> norm < 0.1

##

# repeat with more difficult starting point

initVariable!(fg, :ωβ, [zeros(2) for _ in 1:100])

pts_ = approxConv(fg, :x0x1ωβf1, :ωβ)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [0.7;-0.3]) < 0.1


@warn "n-ary DERelative test on :ωβ requires issue #1010 to be resolved first before being reintroduced."
# ## do a complete solve (must first resolve #1010)

# solveTree!(fg);

# ## Solve quality might not yet be good enough for this particular test case

# @test calcMeanMaxSuggested(fg, :ωβ, :default).suggested - [0.7;-0.3] |> norm < 0.2
# for sym in setdiff(ls(tfg), [:ωβ])
#   @test calcMeanMaxSuggested(fg, sym, :default).suggested - sl(getVariable(fg, sym) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.2
# end


##

end





@error "DERelative not tested for `multihypo=` case yet, see issue #1025"




#