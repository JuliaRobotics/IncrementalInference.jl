
struct ProductNumbers <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
getSample(s::ProductNumbers, N::Int=1) = (reshape(rand(s.z,N),1,:), )

function (s::ProductNumbers)(res::Array{Float64},
      userdata::FactorMetadata,
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      Y::Array{Float64,2},
      XY::Array{Float64,2}  )
  #
  res[1] = XY[1,idx] - X[1,idx]*Y[1,idx] + meas[1][1,idx]
  nothing
end



struct AreEqual <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
function getSample(s::AreEqual, N::Int=1)
  return (reshape(rand(s.z,N),1,:), )
end

function (s::AreEqual)(res::Array{Float64},
                       userdata::FactorMetadata,
                       idx::Int,
                       meas::Tuple{Array{Float64,2}},
                       X::Array{Float64,2},
                       Y::Array{Float64,2}  )
  #
  res[1] = X[1,idx]-Y[1,idx] + meas[1][1,idx]
  nothing
end



struct Square <: IncrementalInference.FunctorPairwise
  z::Distributions.Normal
end
getSample(s::Square, N::Int=1) = (reshape(rand(s.z,N),1,:), )

function (s::Square)(res::Array{Float64},
      userdata::FactorMetadata,
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      X::Array{Float64,2},
      XX::Array{Float64,2}  )
  #
  res[1] = XX[1,idx] - X[1,idx]*X[1,idx] + meas[1][1,idx]
  nothing
end


mutable struct NumbersPrior <: IncrementalInference.FunctorSingleton
  z::BallTreeDensity
end
getSample(s::NumbersPrior, N::Int=1) = (reshape(rand(s.z,N),1,:), )
