
"""
    $SIGNATURES

Return a random sample as a tangent vector from a belief represented by coordinates on a manifold at point p.

Notes

"""
function sampleTangent end

# Sampling MKD
function sampleTangent(M::AbstractDecoratorManifold, x::ManifoldKernelDensity, p = mean(x))
  # get legacy matrix of coordinates and selected labels
  #TODO make sure that when `sample` is replaced in MKD, coordinates is a vector
  coords, lbls = sample(x.belief, 1)
  X = hat(x.manifold, p, coords[:])
  return X
end

function sampleTangent(x::ManifoldKernelDensity, p = mean(x))
  return sampleTangent(x.manifold, x, p)
end

# Sampling Distributions
# assumes M is a group and will break for Riemannian, but leaving that enhancement as TODO
function sampleTangent(
  M::AbstractManifold,
  z::Distribution,
  p = getPointIdentity(M),
  basis::AbstractBasis = DefaultOrthogonalBasis()
)
  return get_vector(M, p, rand(z), basis)
end

function sampleTangent(
  M::AbstractDecoratorManifold,
  z::Distribution,
  p = getPointIdentity(M),
)
  return hat(M, p, SVector{length(z)}(rand(z))) #TODO make sure all Distribution has length, 
                                                # if this errors maybe fall back no next line
  # return convert(typeof(p), hat(M, p, rand(z, 1)[:])) #TODO find something better than (z,1)[:]
end

"""
    $SIGNATURES

Return a random sample point on a manifold from a belief represented by coordinates at point p.

Notes

"""
function samplePoint(
  M::AbstractManifold,
  sbelief,
  p,
  basis::AbstractBasis,
  retraction_method::AbstractRetractionMethod = ExponentialRetraction(),
)
  X = sampleTangent(M, sbelief, p, basis)
  return retract(M, p, X, retraction_method)
end
function samplePoint(
  M::AbstractDecoratorManifold,
  sbelief,
  p = getPointIdentity(M),
  retraction_method::AbstractRetractionMethod = ExponentialRetraction(),
)
  X = sampleTangent(M, sbelief, p)
  return retract(M, p, X, retraction_method)
end

function samplePoint(
  M::AbstractDecoratorManifold,
  sbelief::ManifoldKernelDensity,
  # p = identity_element(M, mean(sbelief)), # 8.671254 seconds (82.64 M allocations: 3.668 GiB, 7.50% gc time)
  p = getPointIdentity(M), #6.713209 seconds (66.42 M allocations: 3.141 GiB, 7.52% gc time)
  retraction_method::AbstractRetractionMethod = ExponentialRetraction(),
)
  X = sampleTangent(M, sbelief, p)
  return retract(M, p, X, retraction_method)
end

function samplePoint(x::ManifoldKernelDensity, p = mean(x))
  return samplePoint(x.manifold, x, p)
end

# FIXME: rather use manifolds
function samplePoint(distr::SamplableBelief)
  Base.depwarn(
    "samplePoint(distr::SamplableBelief) should be replaced by samplePoint(M<:AbstractManifold, distr::SamplableBelief, ...)",
    :samplePoint,
  )
  return rand(distr, 1)
end

## default getSample
"""
    $SIGNATURES

Sample the factor in `CalcFactor`. A default `getSample` method is provided that should cover most use cases, 
if more advanced sampling is required, the `getSample` function should be extended.

The default behavior for `getSample` is as follows:
- The `SamplableBelief`` shall be in the field `Z` and that shall be enough to fully define the factor, i.e. `Z<:SamplableBelief` should be the only field.
- Sampling on `<:AbstractManifoldMinimize` factors defined on Group Manifolds: 
  - `getSample` normally returns a tangent vector at the identity element, however it should just match the custom factor definition.
- Sampling on prior (`<:AbstractPrior`) factors : 
  - `getSample` must return a point on the manifold that matches the point representation of the variable.

Notes
- Users should overload this method should their factor not only use field `Z` for the `SamplableBelief`.
- See the Custom Factors section in the Caesar.jl documentation for more examples and details.
- Also see issue https://github.com/JuliaRobotics/IncrementalInference.jl/issues/1441

See also: [`getMeasurementParametric`](@ref)
"""
function getSample end

function getSample(cf::CalcFactor{<:AbstractPrior})
  M = getManifold(cf)
  if hasfield(typeof(cf.factor), :Z)
    X = samplePoint(M, cf.factor.Z)
  else
    error(
      """Factor $(typeof(cf.factor)) does not have a field `Z`, to use the default `getSample` method, use `Z` for the measurement. 
          Alternatively, provide a `getSample` method. See IIF issue #1441 and Custom Factors in the Caesar documentation.""",
    )
  end
  return X
end

function getSample(cf::CalcFactor{<:AbstractRelative})
  M = getManifold(cf)
  if hasfield(typeof(cf.factor), :Z)
    X = sampleTangent(M, cf.factor.Z)
  else
    error(
      """Factor $(typeof(cf.factor)) does not have a field `Z`, to use the default `getSample` method, use `Z` for the measurement. 
          Alternatively, provide a `getSample` method. See IIF issue #1441 and Custom Factors in the Caesar documentation.""",
    )
  end
  return X
end
