module InterpolationsExt

@info "IncrementalInference.jl is loading extensions related to Interpolations.jl."

using Interpolations
using Statistics
using DocStringExtensions
using TensorCast
using Manifolds
using ApproxManifoldProducts
import ApproxManifoldProducts: sample
const AMP = ApproxManifoldProducts

import IncrementalInference: getManifold
import IncrementalInference: HeatmapGridDensity, PackedHeatmapGridDensity
import IncrementalInference: LevelSetGridNormal, PackedLevelSetGridNormal

export HeatmapGridDensity, PackedHeatmapGridDensity
export LevelSetGridNormal, PackedLevelSetGridNormal
export sampleHeatmap


include("HeatmapSampler.jl")

end # module