# test null hypothesis cases

using Base: Test
using IncrementalInference, Distributions
# going to introduce two new constraint types
import IncrementalInference: getSample


type DevelopPriorNH <: IncrementalInference.FunctorSingletonNH
  x::Distribution
  nullhypothesis::Distributions.Categorical
end
getSample(dpl::DevelopPriorNH, N::Int=1) = (rand(dpl.x, N)', )



N  = 100
fg = emptyFactorGraph()

v1 = addNode!(fg,:x1,ones(1,N),N=N)

pr = DevelopPriorNH(Normal(10.0,1.0), Categorical([0.5;0.5]))
f1 = addFactor!(fg,[:x1],pr)

# Juno.breakpoint("/home/dehann/.julia/v0.5/IncrementalInference/src/ApproxConv.jl",121)

pts = evalFactor2(fg, f1, v1.index, N=N)

@test sum(abs.(pts - 1.0) .< 5) > 30
@test sum(abs.(pts - 10.0) .< 5) > 30






#
