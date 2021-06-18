using Manopt, Manifolds
using LinearAlgebra 
using StaticArrays

## The factors
struct MeasurementOnTangent end

function PointPoint_distance(M, m, p, q)
    q̂ = compose(M, p, m)
    return distance(M, q, q̂)
end

function PointPoint_distance(M, X, p, q, ::MeasurementOnTangent)
    q̂ = compose(M, p, exp(M, identity(M, X), X)) #FIXME
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

# Pose2Pose2(m, p, q) = PointPoint_distance(getManifold(Pose2), m, p, q)
Pose2Pose2(m, p, q) = PointPoint_distance(SpecialEuclidean(2), m, p, q)
Pose3Pose3(m, p, q) = PointPoint_distance(SpecialEuclidean(3), m, p, q)
PriorPose2(m, p) = Prior_distance(SpecialEuclidean(2), m, p)
PriorPose3(m, p) = Prior_distance(SpecialEuclidean(3), m, p)


## Testing the Factors

ϵSE2 = ProductRepr([0.,0] , [1. 0; 0 1])
ϵSE3 = ProductRepr([0.,0,0] , [1. 0 0; 0 1 0; 0 0 1])

M = SpecialEuclidean(2)

X = hat(M, ϵSE2, [1, 0, pi/4])
p = compose(M, ϵSE2, exp(M, ϵSE2, X))

X = hat(M, ϵSE2, [1, 0, pi/8])
q = compose(M, p, exp(M, ϵSE2, X))

X = hat(M, ϵSE2, [1-1e-3, 0+1e-3, pi/8+1e-3])
m = compose(M, ϵSE2, exp(M, ϵSE2, X))

Pose2Pose2(m, p, q)

PointPoint_distance(M, X, p, q, MeasurementOnTangent())

PriorPose2(m, q)
PriorPose2(q, q)


## and with matrices
using LinearAlgebra

ϵSE2 = Matrix(1.0I(3))
ϵSE3 = Matrix(1.0I(4))

M = SpecialEuclidean(2)

X = hat(M, ϵSE2, [1, 0, pi/4])
# exp does not work with affine matrix
p = compose(M, ϵSE2, exp(M, ϵSE2, X))


##

# 2 variables x1 and x2 of type SE2 one factor and one prior on x1

X = hat(M, ϵSE2, [1, 0, pi/4])
priorx1 = compose(M, ϵSE2, exp(M, ϵSE2, X))

measX = hat(M, ϵSE2, [1, 0, pi/4])
measx1x2 = compose(M, ϵSE2, exp(M, ϵSE2, measX))


function F_prior(M,x) 
    return PriorPose2(priorx1, x)
end

function gradF_prior(M, x)
    gradient(M,(x)->F_prior(M,x), x)
end

x = deepcopy(ϵSE2)

F_prior(M, x)
gradF_prior(M, x)

# does something but wrong
xMean = gradient_descent(M, F_prior, gradF_prior, x)

F_prior(M, xMean)

stopping_criterion = (StopWhenAny(StopAfterIteration(1000),StopWhenGradientNormLess(10.0^-8))) 
o = gradient_descent(M, F_prior, gradF_prior, x; stopping_criterion, return_options=true)
o.stop
o.x
#also errors
n = 100
x0 = [Manopt.random_point(M) for i=1:n] 
particle_swarm(M, F_prior; n, x0)


##

MP = ProductManifold(M,M)
function F(MP,x) 
    
    return PriorPose2(priorx1, x.parts[1]) +  Pose2Pose2(measx1x2, x.parts[1], x.parts[2])
end

function gradF(MP, x)
    gradient(M,(x)->F(M,x), x)
end

xp = ProductRepr(deepcopy(ϵSE2), deepcopy(ϵSE2))
F(MP, xp)

# does not work
gradF(MP, xp)



xMean = gradient_descent(MP, F, gradF, xp)


# does something but wrong
NelderMead(MP, F)

#errors 
NelderMead(MP, F, xp)



## ======================================================================================
## with Optim
## ======================================================================================

using Optim

# from IIF

struct ManifoldsVector <: Optim.Manifold
  manis::Vector{AbstractManifold}
end

ManifoldsVector() = ManifoldsVector(AbstractManifold[])

Base.getindex(mv::ManifoldsVector, inds...) = getindex(mv.manis, inds...)
Base.setindex!(mv, X, inds...) =  setindex!(mv.manis, X, inds...)
Base.iterate(mv::ManifoldsVector, args...) = iterate(mv.manis, args...)

# function ManifoldsVector(fg::AbstractDFG, varIds::Vector{Symbol})
#   manis = Bool[]
#   for k = varIds
#     push!(manis, getVariableType(fg, k) |> getManifold)
#   end
#   ManifoldsVector(manis)
# end

function Optim.retract!(manis::ManifoldsVector, x)

    for (i,M) = enumerate(manis)
        _x = uncoords(x[:,i])
        project!(M, _x, _x)
        x[:,i] .= coords(_x)
    end
    return x 
end

function Optim.project_tangent!(manis::ManifoldsVector, G, x)

    for (i, M) = enumerate(manis)
        _x = uncoords(x[:, i])
        _G = hat(M, identity(M, _x), G)
        project!(M, _G, _x, _G)
        G[:, i] .= vee(M, identity(M, _x), _G)
    end
    return G
end





function coords(p)
    return [p.parts[1][1], p.parts[1][2], atan(p.parts[2][2,1],p.parts[2][1,1])]
end
# reverse of `coords`
function uncoords(p)
    α = p[3]
    return ProductRepr(([p[1], p[2]]), [cos(α) -sin(α); sin(α) cos(α)])
end


autodiff = :forward
autodiff = :finite
algorithm=Optim.BFGS
algorithm=Optim.NelderMead
algorithmkwargs=() # add manifold to overwrite computed one
options = Optim.Options(allow_f_increases=true, 
                        iterations = 100,
                        time_limit = 100,
                        show_trace = true,
                        show_every = 1,
                        )
##



coords(priorx1 )

function cost(x) 
    # @show "cost" x
    _x = uncoords(x)
    return PriorPose2(priorx1, _x)
end

manifold = ManifoldsVector([SpecialEuclidean(2)])
alg = algorithm(;manifold, algorithmkwargs...)
# alg = algorithm(; algorithmkwargs...)


initValues = collect([0. 0 0]')
tdtotalCost = Optim.TwiceDifferentiable((x)->cost(x), initValues; autodiff)
tdtotalCost = Optim.NonDifferentiable((x)->cost(x), initValues)


result = Optim.optimize(tdtotalCost, initValues, alg, options)

result = Optim.optimize(cost, initValues, alg, options)

rv = Optim.minimizer(result)
uncoords(rv)




H = Optim.hessian!(tdtotalCost, rv)

Σ = pinv(H)

