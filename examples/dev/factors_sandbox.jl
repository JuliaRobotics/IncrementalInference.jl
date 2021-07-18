using Manifolds
using LinearAlgebra 
using StaticArrays

## The factors
struct MeasurementOnTangent end

function PointPoint_distance(M, m, p, q)
    q̂ = compose(M, p, m)
    return distance(M, q, q̂)
end

function grad_PointPoint_distance(M, m, p, q) 
    q̂ = compose(M, p, m)
    return grad_distance(M, q, q̂)
end

function PointPoint_distance(M, X, p, q, ::MeasurementOnTangent)
    q̂ = compose(M, p, exp(M, identity(M, p), X))
    return distance(M, q, q̂)
end

function PointPoint_velocity_distance(M, X, dt, p, q, ::MeasurementOnTangent)
    q̂ = compose(M, p, group_exp(M, X*dt))
    return distance(M, q, q̂)
end

function Prior_distance(M, meas, p)	
#		
    return distance(M, meas, p)
end

# Pose
# Pose2Pose2(m, p, q) = PointPoint_distance(getManifold(Pose2), m, p, q)
Pose2Pose2(m, p, q) = PointPoint_distance(SpecialEuclidean(2), m, p, q)
Pose3Pose3(m, p, q) = PointPoint_distance(SpecialEuclidean(3), m, p, q)
PriorPose2(m, p) = Prior_distance(SpecialEuclidean(2), m, p)
PriorPose3(m, p) = Prior_distance(SpecialEuclidean(3), m, p)

#Point
Point1Point1(m, p, q) = PointPoint_distance(TranslationGroup(1), m, p, q)
Point2Point2(m, p, q) = PointPoint_distance(TranslationGroup(2), m, p, q)
Point3Point3(m, p, q) = PointPoint_distance(TranslationGroup(3), m, p, q)
PriorPoint2(m, p) = Prior_distance(TranslationGroup(2), m, p)
PriorPoint3(m, p) = Prior_distance(TranslationGroup(3), m, p)


## Testing the Factors

ϵSE2 = ProductRepr(SA[0.,0] , SA[1. 0; 0 1])
ϵSE3 = ProductRepr(SA[0.,0,0] , SA[1. 0 0; 0 1 0; 0 0 1])

M = SpecialEuclidean(2)

X = hat(M, ϵSE2, [1, 0, pi/4+0.1])
p = compose(M, ϵSE2, exp(M, ϵSE2, X))

X = hat(M, ϵSE2, [1, 0, pi/8])
q = compose(M, p, exp(M, ϵSE2, X))

X = hat(M, ϵSE2, [1-1e-3, 0+1e-3, pi/8+1e-3])
m = compose(M, ϵSE2, exp(M, ϵSE2, X))

Pose2Pose2(m, p, q)

PointPoint_distance(M, X, p, q, MeasurementOnTangent())

PriorPose2(m, q)
PriorPose2(q, q)



## Testing the Factor


## and with matrices
# using LinearAlgebra

# ϵSE2 = Matrix(1.0I(3))
# ϵSE3 = Matrix(1.0I(4))

# M = SpecialEuclidean(2)

# X = hat(M, ϵSE2, [1, 0, pi/4])
# # exp does not work with affine matrix
# p = compose(M, ϵSE2, exp(M, ϵSE2, X))