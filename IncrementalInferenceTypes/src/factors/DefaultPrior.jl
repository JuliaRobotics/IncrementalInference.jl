"""
$(TYPEDEF)

Default prior on all dimensions of a variable node in the factor graph.  `Prior` is
not recommended when non-Euclidean dimensions are used in variables.
"""

DFG.@defObservationType Prior AbstractPriorObservation obs->IncrementalInferenceTypes.TranslationGroup(getDimension(obs.Z))

# struct Prior{T} <: AbstractPriorObservation
#   Z::T
# end
# DFG.getManifold(pr::Prior) = LieGroups.TranslationGroup(getDimension(pr.Z))

# """
# $(TYPEDEF)

# Serialization type for Prior.
# """
# Base.@kwdef mutable struct PackedPrior <: AbstractPackedObservation
#   Z::PackedBelief
# end

# function DFG.pack(d::Prior)
#   return PackedPrior(DFG.pack(d.Z))
# end

# function DFG.unpack(d::PackedPrior)
#   return Prior(DFG.unpack(d.Z))
# end

#
