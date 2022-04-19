

# NOTE, user variables and manifolds will require the same definitions, TODO perhaps add into `@defVariable`
# unusual definitions, but all they need to do is pack and unpack as one-to-one
# this step actually occurs separate from the actual variables or factors (with their own manifolds) 
# relies on later use of getManifold to give back the same <:AbstractManifold
# NOTE added to DFG.@defVariable
getVariableType(M::Euclidean{Tuple{N}}) where N = ContinuousEuclid(N)
getVariableType(M::TranslationGroup{Tuple{N}}) where N = ContinuousEuclid(N)

# getVariableType(M::RealCircleGroup) = Circular()
# getVariableType(M::Circle) = error("Circle manifold is deprecated use RealCircleGroup, will come back when we generalize to non-group Riemannian")



# Type converters for MKD
Base.convert(::Type{<:SamplableBelief}, ::Type{<:PackedManifoldKernelDensity}) = ManifoldKernelDensity
Base.convert(::Type{<:PackedSamplableBelief}, ::Type{<:ManifoldKernelDensity}) = PackedManifoldKernelDensity


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
  bw = getBW(mkd)[:,1]
  manikde!(mkd.manifold, pts[1:1], mkd._u0; bw, partial=mkd._partial, infoPerCoord=mkd.infoPerCoord)
end


# Data converters for MKD
function packDistribution( mkd::ManifoldKernelDensity )
  #
  pts = getPoints(mkd)

  PackedManifoldKernelDensity(
    "IncrementalInference.PackedManifoldKernelDensity",
    # piggy back on InferenceVariable serialization rather than try serialize anything Manifolds.jl
    DFG.typeModuleName(getVariableType(mkd.manifold)),
    [AMP.makeCoordsFromPoint(mkd.manifold, pt) for pt in pts],
    getBW(mkd.belief)[:,1],
    mkd._partial isa Nothing ? collect(1:manifold_dimension(mkd.manifold)) : mkd._partial ,
    mkd.infoPerCoord
  )
end


function unpackDistribution(dtr::PackedManifoldKernelDensity)
  # find InferenceVariable type from string (anything Manifolds.jl?)
  M = DFG.getTypeFromSerializationModule(dtr.varType) |> getManifold
  vecP = [AMP.makePointFromCoords(M, pt) for pt in dtr.pts]

  partial = length(dtr.partial) == manifold_dimension(M) ? nothing : dtr.partial
  
  return manikde!( M, vecP; bw=dtr.bw, partial, infoPerCoord=dtr.infoPerCoord )
end



function Base.convert(::Type{String}, 
                      mkd::ManifoldKernelDensity )
  #
  packedMKD = packDistribution(mkd)
  JSON2.write(packedMKD)
end


# Use general dispatch
# Base.convert(::Type{<:PackedSamplableBelief}, mkd::ManifoldKernelDensity) = convert(String, mkd)

# make module specific
# good references: 
#  https://discourse.julialang.org/t/converting-string-to-datatype-with-meta-parse/33024/2
#  https://discourse.julialang.org/t/is-there-a-way-to-import-modules-with-a-string/15723/6
function Base.convert(::Type{<:ManifoldKernelDensity}, str::AbstractString)
  dtr = JSON2.read(str, PackedManifoldKernelDensity)
  unpackDistribution(dtr)
end



#
