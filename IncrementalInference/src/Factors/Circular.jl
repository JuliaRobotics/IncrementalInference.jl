# ---------------------------------------------
# CircularCircular
# ---------------------------------------------

function (cf::CalcFactor{<:CircularCircular})(X, p, q)
  M = getManifold(cf)
  return measurement_residual(M, X, p, q)
end

function getSample(cf::CalcFactor{<:CircularCircular})
  # FIXME workaround for issue with manifolds CircularGroup, 
  # return [rand(cf.factor.Z)]
  return sampleTangent(getManifold(cf), cf.factor.Z, getPointIdentity(Circular))
end

function Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{CircularCircular})
  return Manifolds.RealCircleGroup()
end

IIFTypes.CircularCircular(::UniformScaling) = CircularCircular(Normal())

# ---------------------------------------------
# PriorCircular
# ---------------------------------------------

IIFTypes.PriorCircular(::UniformScaling) = PriorCircular(Normal())

function getSample(cf::CalcFactor{<:PriorCircular})
  # FIXME workaround for issue #TBD with manifolds CircularGroup, 
  # JuliaManifolds/Manifolds.jl#415
  # no method similar(::Float64, ::Type{Float64})
  # return samplePoint(cf.manifold, cf.factor.Z, [0.0])
  return samplePoint(cf.manifold, cf.factor.Z, getPointIdentity(Circular))
  # return [Manifolds.sym_rem(rand(cf.factor.Z))]
end

function (cf::CalcFactor{<:PriorCircular})(m, p)
  M = getManifold(cf)
  Xc = vee(LieAlgebra(M), log(M, p, m))
  return Xc
end

function Base.convert(::Type{<:MB.AbstractManifold}, ::InstanceType{PriorCircular})
  return Manifolds.RealCircleGroup()
end

# --------------------------------------------
