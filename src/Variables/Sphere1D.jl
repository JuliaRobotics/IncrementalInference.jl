
export Sphere1

"""
$(TYPEDEF)

Sphere1 is a S1 mechanization of one Circular rotation, with `theta in [-pi,pi)`.
"""
struct Sphere1 <: InferenceVariable
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol}
  Sphere1(;labels::Vector{<:AbstractString}=String[]) = new(1, labels, (:Circular,))
end
