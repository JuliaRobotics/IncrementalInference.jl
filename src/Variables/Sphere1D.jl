
export Sphere1

"""
$(TYPEDEF)

Sphere1 is a S1 mechanization of one Circular rotation, with `theta in [-pi,pi)`.
"""
struct Sphere1 <: InferenceVariable
  dims::Int
  manifolds::Tuple{Symbol}
  Sphere1() = new(1, (:Circular,))
end
