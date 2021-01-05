
struct ProductNumbers <: AbstractRelativeRoots
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



struct AreEqual <: AbstractRelativeRoots
  z::Distributions.Normal
end
function getSample(cf::CalcFactor{<:AreEqual}, N::Int=1)
  return (reshape(rand(cf.factor.z,N),1,:), )
end

function (cf::CalcFactor{<:AreEqual})( res::AbstractVector{<:Real},
                                      meas,
                                      X,
                                      Y  )
  #
  res[1] = X[1]-Y[1] + meas[1]
  nothing
end



struct Square <: AbstractRelativeRoots
  z::Distributions.Normal
end
getSample(cf::CalcFactor{<:Square}, N::Int=1) = (reshape(rand(cf.factor.z,N),1,:), )

function (cf::CalcFactor{<:Square})( res::AbstractVector{<:Real},
                                    meas,
                                    X,
                                    XX  )
  #
  res[1] = XX[1] - X[1]*X[1] + meas[1]
  nothing
end


mutable struct NumbersPrior <: AbstractPrior
  z::BallTreeDensity
end
getSample(cf::CalcFactor{<:NumbersPrior}, N::Int=1) = (reshape(rand(cf.factor.z,N),1,:), )
