
struct ProductNumbers <: AbstractRelativeFactor
  z::Distributions.Normal
end
getSample(s::ProductNumbers, N::Int=1) = (reshape(rand(s.z,N),1,:), )

function (s::ProductNumbers)(res::AbstractArray{<:Real},
      userdata::FactorMetadata,
      idx::Int,
      meas::Tuple{<:AbstractArray{<:Real,2}},
      X::AbstractArray{<:Real,2},
      Y::AbstractArray{<:Real,2},
      XY::AbstractArray{<:Real,2}  )
  #
  res[1] = XY[1,idx] - X[1,idx]*Y[1,idx] + meas[1][1,idx]
  nothing
end



struct AreEqual <: AbstractRelativeFactor
  z::Distributions.Normal
end
function getSample(s::AreEqual, N::Int=1)
  return (reshape(rand(s.z,N),1,:), )
end

function (s::AreEqual)(res::AbstractArray{<:Real},
                       userdata::FactorMetadata,
                       idx::Int,
                       meas::Tuple{<:AbstractArray{<:Real,2}},
                       X::AbstractArray{<:Real,2},
                       Y::AbstractArray{<:Real,2}  )
  #
  res[1] = X[1,idx]-Y[1,idx] + meas[1][1,idx]
  nothing
end



struct Square <: AbstractRelativeFactor
  z::Distributions.Normal
end
getSample(s::Square, N::Int=1) = (reshape(rand(s.z,N),1,:), )

function (s::Square)(res::AbstractArray{<:Real},
      userdata::FactorMetadata,
      idx::Int,
      meas::Tuple{<:AbstractArray{<:Real,2}},
      X::AbstractArray{<:Real,2},
      XX::AbstractArray{<:Real,2}  )
  #
  res[1] = XX[1,idx] - X[1,idx]*X[1,idx] + meas[1][1,idx]
  nothing
end


mutable struct NumbersPrior <: AbstractPrior
  z::BallTreeDensity
end
getSample(s::NumbersPrior, N::Int=1) = (reshape(rand(s.z,N),1,:), )
