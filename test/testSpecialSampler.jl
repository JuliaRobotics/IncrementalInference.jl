# test new sampling API

using Test
using IncrementalInference

# import IncrementalInference: getSample


struct SpecialPrior{T} <: IncrementalInference.FunctorSingleton where T <: Distribution
  z::T
  specialSampler::Function
end
specialSample(s::SpecialPrior, N::Int, fmd::FactorMetadata, Xi) = (reshape(rand(s.z,N),1,:), )
SpecialPrior(z::T) where {T <: SamplableBelief}= SpecialPrior{T}(z, specialSample)

struct SpecialLinearOffset{T} <: IncrementalInference.FunctorPairwise where T <: Distribution
  z::T
  specialSampler::Function
end

specialSample(s::SpecialLinearOffset, N::Int, fmd::FactorMetadata, Xi, Xj) = (reshape(rand(s.z,N),1,:), )

SpecialLinearOffset(z::T) where {T <: SamplableBelief} = SpecialLinearOffset{T}(z, specialSample)

function (s::SpecialLinearOffset)(res::Array{Float64},
                           userdata::FactorMetadata,
                           idx::Int,
                           meas::Tuple{Array{Float64, 2}},
                           X1::Array{Float64,2},
                           X2::Array{Float64,2}  )
  #
  res[1] = meas[1][idx] - (X2[1,idx] - X1[1,idx])
  nothing
end


@testset "test specialSampler functionality..." begin

fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0], SpecialPrior(Normal()))

addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0;:x1], SpecialLinearOffset(Normal(10,1)))

tree, smt, hist = solveTree!(fg)

@test getPPE(fg, :x0).suggested[1] |> abs < 1.0

@test getPPE(fg, :x1).suggested[1] - 10 |> abs < 3.0


end



## specific test for #568



#
