# test null hypothesis cases

using Test
using IncrementalInference
# going to introduce two new constraint types
import IncrementalInference: getSample


mutable struct DevelopPriorNH <: IncrementalInference.FunctorSingletonNH
  x::Distribution
  nullhypothesis::Distributions.Categorical
end
function getSample(dpl::DevelopPriorNH, N::Int=1)
  return (reshape(rand(dpl.x, N),1,N), )
end


@testset "test null hypothesis singletons..." begin

global N  = 100
global fg = initfg()

global v1 = addVariable!(fg, :x1, ContinuousScalar, N=N)

global pr = DevelopPriorNH(Normal(10.0,1.0), Categorical([0.5;0.5]))
global f1 = addFactor!(fg,[:x1],pr, autoinit=true)

# ensureAllInitialized!(fg)
# Juno.breakpoint("/home/dehann/.julia/v0.5/IncrementalInference/src/ApproxConv.jl",121)

global pts = evalFactor2(fg, f1, v1.label, N=N)

@test sum(abs.(pts .- 1.0) .< 5) > 30
@test sum(abs.(pts .- 10.0) .< 5) > 30

end



#
