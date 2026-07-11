
IIFTypes.Prior(::UniformScaling) = Prior(Normal())

# getSample(cf::CalcFactor{<:Prior}) = rand(cf.factor.Z)

# basic default
(s::CalcFactor{<:Prior})(z, x1) = z .- x1
