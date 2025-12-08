
"""
$(TYPEDEF)

Default linear offset between two scalar variables.

```math
X_2 = X_1 + η_Z
```
"""
# @tags struct LinearRelative{T} <: RelativeObservation
#   Z::T & DFG.@packed
# end

#TODO
DFG.@defObservationType LinearRelative RelativeObservation obs->IncrementalInferenceTypes.TranslationGroup(getDimension(obs.Z))

# DFG.getManifold(obs::LinearRelative) = LieGroups.TranslationGroup(getDimension(obs.Z))

# """
# $(TYPEDEF)
# Serialization type for `LinearRelative` binary factor.
# """
# Base.@kwdef mutable struct PackedLinearRelative <: AbstractPackedObservation
#   Z::PackedBelief
# end

# function DFG.pack(d::LinearRelative)
#   return PackedLinearRelative(DFG.pack(d.Z))
# end

# function DFG.unpack(d::PackedLinearRelative)
#   return LinearRelative(DFG.unpack(d.Z))
# end

#
