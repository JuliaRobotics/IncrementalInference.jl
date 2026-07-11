function Base.convert(
  ::Type{<:ManifoldsBase.AbstractManifold},
  ::InstanceType{Position{N}},
) where {N}
  return TranslationGroup(N)
end