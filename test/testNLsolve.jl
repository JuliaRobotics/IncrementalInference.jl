using Test
using NLsolve
using Distributions
# using Statistics

using IncrementalInference

import IncrementalInference: getSample


##

@testset "Ensure lambda's work with NLsolve" begin

function f!(F, x, dummy)
    F[1] = (x[1]+3)*(x[2]^3-7)+18
    F[2] = sin(x[2]*exp(x[1])-1)
end

data = randn(2,3)
r = nlsolve( (res, x) -> f!(res, x, data), zeros(2), inplace=true )

@test norm(r.zero - [0, 1.0]) < 1e-10

end



## abstract functional defintion

mutable struct OneDimensionTest{T} <: AbstractRelativeRoots
  Dx::T
  OneDimensionTest{T}(a::T) where T = new(a)
end
OneDimensionTest(a::T) where T = OneDimensionTest{T}(a)

getSample(od::OneDimensionTest, N::Int=1) = (reshape(rand(od.Dx,N),1,:),)

# standardized parameter list used by IIF
function (Dp::OneDimensionTest)(res::AbstractArray{<:Real},
                                userdata::FactorMetadata,
                                idx::Int,
                                meas::Tuple{<:AbstractArray{<:Real,2}},
                                p1::AbstractArray{<:Real,2},
                                p2::AbstractArray{<:Real,2} )
  #
  res[1] = meas[1][1,idx] - (p2[1,idx] - p1[1,idx])
  nothing
end

##


@testset "minimalistic pattern on how NLsolve is used in IIF." begin

##

fg = initfg()
X0 = addVariable!(fg, :x0, ContinuousScalar)
X1 = addVariable!(fg, :x1, ContinuousScalar)

N = 10
p1 = rand(1,N)
p2 = rand(1,N)
VARS = Array{Array{Float64,2},1}()
push!(VARS,p1)
push!(VARS,p2)

odo = OneDimensionTest(Normal(10.0,1.0))

addFactor!(fg, [:x0;:x1], odo, graphinit=false)


##

meas = (reshape(rand(odo.Dx,N),1,N),)
fmd = IIF._defaultFactorMetadata([X0;X1], arrRef=VARS)
ccw = CommonConvWrapper(odo, VARS[1], 1, VARS, fmd, measurement=meas, varidx=2)

##

# Find definition of (::ccw)(res, x) here:
nlsolve(  ccw, [0.0], inplace=true )

##

end




#
