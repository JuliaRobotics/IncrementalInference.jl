
"""
$(TYPEDEF)

Default linear offset between two scalar variables.

```math
X_2 = X_1 + η_Z
```
"""
struct LinearRelative{T} <: RelativeObservation
  Z::T
end

DFG.getManifold(obs::LinearRelative) = LieGroups.TranslationGroup(getDimension(obs.Z))

"""
$(TYPEDEF)
Serialization type for `LinearRelative` binary factor.
"""
Base.@kwdef mutable struct PackedLinearRelative <: AbstractPackedFactor
  Z::PackedBelief
end

function DFG.pack(d::LinearRelative)
  return PackedLinearRelative(DFG.packDistribution(d.Z))
end

function DFG.unpack(d::PackedLinearRelative)
  return LinearRelative(DFG.unpackDistribution(d.Z))
end

#
