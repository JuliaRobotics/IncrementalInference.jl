
export Sphere1

"""
$(TYPEDEF)

Sphere1 is a S1 mechanization of one Circular rotation, with `theta in [-pi,pi)`.
"""
struct Sphere1 <: InferenceVariable end
getDimension(::Sphere1) = 1
getManifolds(::Sphere1) = (:Circular,)
