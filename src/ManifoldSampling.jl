export sampleTangent
export samplePoint
"""
    $SIGNATURES

Return a random sample as a tangent vector from a belief represented by coordinates on a manifold at point p.

Notes

"""
function sampleTangent end

# Sampling MKD
function sampleTangent(M::AbstractGroupManifold, x::ManifoldKernelDensity, p=mean(x))
    # get legacy matrix of coordinates and selected labels
    coords, lbls = sample(x.belief,1)
    X = hat(x.manifold, p, coords)
    return X
end

function sampleTangent(x::ManifoldKernelDensity, p=mean(x)) 
  return sampleTangent(x.manifold, x, p)
end

# Sampling Distributions
function sampleTangent(M::AbstractManifold, z::Distribution, p, basis::AbstractBasis) 
  return get_vector(M, p, rand(z), basis)
end

function sampleTangent(M::AbstractGroupManifold, z::Distribution, p=identity_element(M)) 
  return hat(M, p, rand(z,1)[:]) #TODO find something better than (z,1)[:]
end


"""
    $SIGNATURES

Return a random sample point on a manifold from a belief represented by coordinates at point p.

Notes

"""
function samplePoint(M::AbstractManifold, sbelief, p, basis::AbstractBasis, retraction_method::AbstractRetractionMethod=ExponentialRetraction())
  X = sampleTangent(M, sbelief, p, basis)
  return retract(M, p, X, retraction_method)
end
function samplePoint(M::AbstractGroupManifold, sbelief, p=identity_element(M), retraction_method::AbstractRetractionMethod=ExponentialRetraction())
  X = sampleTangent(M, sbelief, p)
  return retract(M, p, X, retraction_method)
end

function samplePoint(M::AbstractGroupManifold, sbelief::ManifoldKernelDensity, p=identity_element(M, mean(sbelief)), retraction_method::AbstractRetractionMethod=ExponentialRetraction())
  X = sampleTangent(M, sbelief, p)
  return retract(M, p, X, retraction_method)
end

function samplePoint(x::ManifoldKernelDensity, p=mean(x)) 
  return samplePoint(x.manifold, x, p)
end

# FIXME: rather use manifolds
function samplePoint(distr::SamplableBelief)
  Base.depwarn("samplePoint(distr::SamplableBelief) should be replaced by samplePoint(M<:AbstractManifold, distr::SamplableBelief, ...)", :samplePoint)
  rand(distr,1)
end

## default getSample
# getSample(cf::CalcFactor{<:AbstractPrior}, N::Int=1) = ([samplePoint(getManifold(cf.factor), cf.factor.Z, ) for _=1:N], )
# getSample(cf::CalcFactor{<:AbstractRelative}, N::Int=1) = ([sampleTangent(getManifold(cf.factor), cf.factor.Z) for _=1:N], )