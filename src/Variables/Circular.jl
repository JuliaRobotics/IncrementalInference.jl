
export Circular

"""
$(TYPEDEF)

Sphere1 is a S1 mechanization of one Circular rotation, with `theta in [-pi,pi)`.
"""
@defVariable Circular 1 (:Circular,)


Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{Circular}) = Manifolds.Circle{ℝ}

#