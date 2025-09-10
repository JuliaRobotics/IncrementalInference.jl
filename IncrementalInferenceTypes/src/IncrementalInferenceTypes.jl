module IncrementalInferenceTypes

using DistributedFactorGraphs
using DocStringExtensions
import Manifolds
using LieGroups
using LieGroups: TranslationGroup, ℝ
using Distributions
using StaticArrays
import StructTypes

using Dates: now
using Distributed: nprocs
# using RecursiveArrayTools

# TODO using all interal functions of DFG as a transition step, remove true
DFG.@usingDFG true
using DistributedFactorGraphs: getManifold

# export variable types
export 
    Position,
    Position1,
    Position2,
    Position3,
    Position4,
    ContinuousScalar,
    ContinuousEuclid,
    Circular

#export factor types
export 
    Prior,
    PackedPrior,
    # LinearRelative,
    # PackedLinearRelative,
    CircularCircular,
    PriorCircular,
    PackedCircularCircular,
    PackedPriorCircular

# export packed distributions
export
    PackedCategorical,
    PackedUniform,
    PackedNormal,
    PackedZeroMeanDiagNormal,
    PackedZeroMeanFullNormal,
    PackedDiagNormal,
    PackedFullNormal,
    PackedRayleigh

export
    SolverParams

const IIFTypes = IncrementalInferenceTypes
export IIFTypes

# Variable Definitions
include("variables/DefaultVariableTypes.jl")

# Factor Definitions
include("factors/DefaultPrior.jl")
#FIXME maybe upgrade linear relative to this
# include("factors/LinearRelative.jl")
include("factors/Circular.jl")

# Distribution Serialization
include("serialization/entities/SerializingDistributions.jl")
include("serialization/services/SerializingDistributions.jl")

# solver params
include("solverparams/SolverParams.jl")

end # module IncrementalInferenceTypes
