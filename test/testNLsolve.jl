using Test
using NLsolve
using Distributions
# using Statistics

using IncrementalInference


@testset "Ensure lambda's work with NLsolve" begin

function f!(F, x, dummy)
    F[1] = (x[1]+3)*(x[2]^3-7)+18
    F[2] = sin(x[2]*exp(x[1])-1)
end

data = randn(2,3)
r = nlsolve( (res, x) -> f!(res, x, data), zeros(2), inplace=true )

@test norm(r.zero - [0, 1.0]) < 1e-10

end



# abstract functional defintion
mutable struct OneDimensionTest{T} <: FunctorPairwise
  Dx::T
  OneDimensionTest{T}(a::T) where T = new(a)
end
OneDimensionTest(a::T) where T = OneDimensionTest{T}(a)

# standardized parameter list used by IIF
function (Dp::OneDimensionTest)(res::Array{Float64},
                                userdata::FactorMetadata,
                                idx::Int,
                                meas::Tuple{Array{Float64,2}},
                                p1::Array{Float64,2},
                                p2::Array{Float64,2} )
  #
  res[1] = meas[1][1,idx] - (p2[1,idx] - p1[1,idx])
  nothing
end



@testset "minimalistic pattern on how NLsolve is used in IIF." begin

N = 10
global p1 = rand(1,N)
global p2 = rand(1,N)
global VARS = Array{Array{Float64,2},1}()
push!(VARS,p1)
push!(VARS,p2)

global odo = OneDimensionTest(Normal(10.0,1.0))
global meas = (reshape(rand(odo.Dx,N),1,N),)
global ccw = CommonConvWrapper(odo, VARS[1], 1, VARS, measurement=meas, varidx=2)

# Find definition of (::ccw)(res, x) here:
nlsolve(  ccw, [0.0], inplace=true )

end




#
