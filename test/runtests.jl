using IncrementalInference, KernelDensityEstimate
using Base.Test

include("fourdoortest.jl")
@test true


# data structure conversion tests for protobuffing
println("Testing conversion to packed data structure and back")
pd = convert(PackedVariableNodeData,fg.v[1].attributes["data"])
unpckd = convert(VariableNodeData, pd)

@test compare(fg.v[1].attributes["data"], unpckd)
