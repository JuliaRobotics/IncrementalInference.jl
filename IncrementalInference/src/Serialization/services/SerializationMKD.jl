
# NOTE, user variables and manifolds will require the same definitions, TODO perhaps add into `@defVariable`
# unusual definitions, but all they need to do is pack and unpack as one-to-one
# this step actually occurs separate from the actual variables or factors (with their own manifolds) 
# relies on later use of getManifold to give back the same <:AbstractManifold
# NOTE added to DFG.@defVariable
getVariableType(M::Euclidean{TypeParameter{Tuple{N}}}) where {N} = ContinuousEuclid(N)
getVariableType(M::LieGroups.TranslationGroup{ℝ, TypeParameter{Tuple{N}}}) where {N} = ContinuousEuclid(N)

# getVariableType(M::RealCircleGroup) = Circular()
# getVariableType(M::Circle) = error("Circle manifold is deprecated use RealCircleGroup, will come back when we generalize to non-group Riemannian")

# Type converters for MKD
function Base.convert(::Type{<:SamplableBelief}, ::Type{<:PackedManifoldKernelDensity})
  return ManifoldKernelDensity
end
function Base.convert(::Type{<:PackedBelief}, ::Type{<:ManifoldKernelDensity})
  return PackedManifoldKernelDensity
end

"""
    $SIGNATURES

Parching refers to drying/hollowing out the object to reduce its memory size to the bare minimum.

Notes
- Likely to be used in combination with [stashing](@ref section_stash_unstash) where large data blobs are independently stored.
- For example, a point cloud stored as a MKD is very large and likely to duplicate the already stored point cloud object values,
  - When storing the MKD object, it might make sense to parch the MKD first and just persisting the context of the data.  
  - Reconstituting a full MKD object would then require the inverse, where the parched shell is sized and filled from a separate large data blob.
"""
function parchDistribution(mkd::ManifoldKernelDensity)
  pts = getPoints(mkd)
  bw = getBW(mkd)[:, 1]
  return manikde!(
    mkd.manifold,
    pts[1:1],
    mkd._u0;
    bw,
    partial = mkd._partial,
    infoPerCoord = mkd.infoPerCoord,
  )
end

# Data converters for MKD
function DFG.packDistribution(mkd::ManifoldKernelDensity)
  #
  pts = getPoints(mkd)

  return PackedManifoldKernelDensity(
    "IncrementalInference.PackedManifoldKernelDensity",
    # piggy back on StateType serialization rather than try serialize anything Manifolds.jl
    DFG.typeModuleName(getVariableType(mkd.manifold)),
    [AMP.makeCoordsFromPoint(mkd.manifold, pt) for pt in pts],
    getBW(mkd.belief)[:, 1],
    mkd._partial isa Nothing ? collect(1:manifold_dimension(mkd.manifold)) : mkd._partial,
    mkd.infoPerCoord,
  )
end

function DFG.unpackDistribution(dtr::PackedManifoldKernelDensity)
  # find StateType type from string (anything Manifolds.jl?)
  M = DFG.getTypeFromSerializationModule(dtr.varType) |> getManifold
  vecP = [AMP.makePointFromCoords(M, pt) for pt in dtr.pts]
  bw = length(dtr.bw) === 0 ? nothing : dtr.bw
  partial = if length(dtr.partial) == manifold_dimension(M) || length(dtr.partial) === 0
    nothing
  else
    dtr.partial
  end

  return manikde!(M, vecP; bw, partial, infoPerCoord = dtr.infoPerCoord)
end

function Base.convert(::Type{String}, mkd::ManifoldKernelDensity)
  #
  packedMKD = packDistribution(mkd)
  return JSON3.write(packedMKD)
end

# Use general dispatch
# Base.convert(::Type{<:PackedBelief}, mkd::ManifoldKernelDensity) = convert(String, mkd)

# make module specific
# good references: 
#  https://discourse.julialang.org/t/converting-string-to-datatype-with-meta-parse/33024/2
#  https://discourse.julialang.org/t/is-there-a-way-to-import-modules-with-a-string/15723/6
function Base.convert(::Type{<:ManifoldKernelDensity}, str::AbstractString)
  dtr = JSON3.read(str, PackedManifoldKernelDensity)
  return unpackDistribution(dtr)
end

#
