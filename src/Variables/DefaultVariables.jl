
## Euclid 1

"""
    $TYPEDEF

Continuous Euclidean variable of dimension `N` representing a Position in cartesian space.
"""
struct Position{N} <: InferenceVariable end

Position(N::Int) = Position{N}()

# not sure if these overloads are necessary since DFG 775?
DFG.getManifold(::InstanceType{Position{N}}) where {N} = TranslationGroup(N)
function DFG.getDimension(val::InstanceType{Position{N}}) where {N}
  return manifold_dimension(getManifold(val))
end
DFG.getPointType(::Type{Position{N}}) where {N} = Vector{Float64}
DFG.getPointIdentity(M_::Type{Position{N}}) where {N} = zeros(N) # identity_element(getManifold(M_), zeros(N)) 

function Base.convert(
  ::Type{<:ManifoldsBase.AbstractManifold},
  ::InstanceType{Position{N}},
) where {N}
  return TranslationGroup(N)
end

#

"""
$(TYPEDEF)

Most basic continuous scalar variable in a `::DFG.AbstractDFG` object.

Alias of `Position{1}`
"""
const ContinuousScalar = Position{1}
const ContinuousEuclid{N} = Position{N}

const Position1 = Position{1}
const Position2 = Position{2}
const Position3 = Position{3}
const Position4 = Position{4}

## Circular

"""
$(TYPEDEF)

Circular is a `Manifolds.Circle{ℝ}` mechanization of one rotation, with `theta in [-pi,pi)`.
"""
@defVariable Circular RealCircleGroup() [0.0;]
#TODO This is an example of what we want working, possible issue upstream in Manifolds.jl
# @defVariable Circular RealCircleGroup() Scalar(0.0)

#
