
## ================================================================================================================================
# Serialization
## ================================================================================================================================

# abstract type JSONManifoldKernelDensity end

import DistributedFactorGraphs: getVariableType

export PackedManifoldKernelDensity

# NOTE, user variables and manifolds will require the same definitions, TODO perhaps add into `@defVariable`
# unusual definitions, but all they need to do is pack and unpack as one-to-one
# this step actually occurs separate from the actual variables or factors (with their own manifolds) 
# relies on later use of getManifold to give back the same <:AbstractManifold
# NOTE added to DFG.@defVariable
getVariableType(M::Euclidean{Tuple{N}}) where N = ContinuousEuclid(N)
getVariableType(M::TranslationGroup{Tuple{N}}) where N = ContinuousEuclid(N)

# getVariableType(M::RealCircleGroup) = Circular()
# getVariableType(M::Circle) = error("Circle manifold is deprecated use RealCircleGroup, will come back when we generalize to non-group Riemannian")



##

struct PackedManifoldKernelDensity <: PackedSamplableBelief
  _type::String
  varType::String
  pts::Vector{Vector{Float64}}
  bw::Vector{Float64}
  partial::Vector{Int}
  infoPerCoord::Vector{Float64}
end

# Type converters for MKD
Base.convert(::Type{<:SamplableBelief}, ::Type{<:PackedManifoldKernelDensity}) = ManifoldKernelDensity
Base.convert(::Type{<:PackedSamplableBelief}, ::Type{<:ManifoldKernelDensity}) = PackedManifoldKernelDensity

# Data converters for MKD
function Base.convert(::Type{<:AbstractString}, 
                      mkd::ManifoldKernelDensity )
  #
  pts = getPoints(mkd)

  packedMKD = PackedManifoldKernelDensity(
    "PackedManifoldKernelDensity",
    # piggy back on InferenceVariable serialization rather than try serialize anything Manifolds.jl
    DFG.typeModuleName(getVariableType(mkd.manifold)),
    [AMP.makeCoordsFromPoint(mkd.manifold, pt) for pt in pts],
    getBW(mkd.belief)[:,1],
    mkd._partial isa Nothing ? collect(1:manifold_dimension(mkd.manifold)) : mkd._partial ,
    mkd.infoPerCoord
  )
  
  JSON2.write(packedMKD)
end

Base.convert(::Type{<:PackedSamplableBelief}, mkd::ManifoldKernelDensity) = convert(String, mkd)


# make module specific
# good references: 
#  https://discourse.julialang.org/t/converting-string-to-datatype-with-meta-parse/33024/2
#  https://discourse.julialang.org/t/is-there-a-way-to-import-modules-with-a-string/15723/6
function Base.convert(::Type{<:ManifoldKernelDensity}, str::AbstractString)
  data = JSON2.read(str, PackedManifoldKernelDensity)
  
  # piggy back on serialization of InferenceVariable rather than try serialize anything Manifolds.jl
  M = DFG.getTypeFromSerializationModule(data.varType) |> getManifold
  vecP = [AMP.makePointFromCoords(M, pt) for pt in data.pts]
  partial = length(data.partial) == manifold_dimension(M) ? nothing : data.partial
  
  return manikde!( M, vecP, bw=data.bw, partial=partial, infoPerCoord=data.infoPerCoord )
end



#
