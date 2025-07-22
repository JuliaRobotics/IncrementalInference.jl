#FIXME This is discouraged in the julia style guide, rather standardize to instance or type
const InstanceType{T} = Union{Type{<:T}, <:T}

## Euclid 1

"""
    $TYPEDEF

Continuous Euclidean variable of dimension `N` representing a Position in cartesian space.
"""
struct Position{N} <: StateType{N} end

Position(N::Int) = Position{N}()

# not sure if these overloads are necessary since DFG 775?
DFG.getManifold(::InstanceType{Position{N}}) where {N} = LieGroups.TranslationGroup(N)
function DFG.getDimension(val::InstanceType{Position{N}}) where {N}
  return manifold_dimension(getManifold(val))
end
DFG.getPointType(::Type{Position{N}}) where {N} = SVector{N, Float64}
DFG.getPointIdentity(M_::Type{Position{N}}) where {N} = @SVector(zeros(N)) # identity_element(getManifold(M_), zeros(N)) 


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

#TODO maybe just use @defVariable for all Position types?
# @defVariable Position1 TranslationGroup(1) @SVector(zeros(1))
# @defVariable Position2 TranslationGroup(2) @SVector(zeros(2))
# @defVariable Position3 TranslationGroup(3) @SVector(zeros(3))
# @defVariable Position4 TranslationGroup(4) @SVector(zeros(4))

## Circular

"""
$(TYPEDEF)

Circular is a `Manifolds.Circle{ℝ}` mechanization of one rotation, with `theta in [-pi,pi)`.
"""
# @defVariable Circular CircleGroup(ℝ) [0.0;]
# @defVariable(Circular, ValidationLieGroup(LieGroups.CircleGroup()), Scalar(1.0 + 0.0im))
# @defVariable(Circular, LieGroups.CircleGroup(), Scalar(1.0 + 0.0im))
@defVariable(Circular, LieGroups.CircleGroup(), fill(1.0 + 0.0im))
# @defVariable Circular LieGroups.CircleGroup(ℝ) 0.0
#TODO This is an example of what we want working, possible issue upstream in Manifolds.jl
# @defVariable Circular CircleGroup(ℝ) Scalar(0.0)

#
