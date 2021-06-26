include("factors_sandbox.jl")

using Optim

# TODO 2 variables x1 and x2 of type SE2 one factor and one prior on x1

# X = hat(M, ϵSE2, [1.1, 0.1, pi+0.1])
X = hat(M, ϵSE2, [1.1, 0.1, pi + 0.05])
priorx1 = compose(M, ϵSE2, exp(M, ϵSE2, X))

# X = hat(M, ϵSE2, [1.1, 0.1, pi-0.1])
X = hat(M, ϵSE2, [0.9, -0.1, pi + 0.15])
priorx2 = compose(M, ϵSE2, exp(M, ϵSE2, X))

measX = hat(M, ϵSE2, [1, 0, pi/4])
measx1x2 = compose(M, ϵSE2, exp(M, ϵSE2, measX))

X2 = hat(M, ϵSE2, [0, 0, pi + pi/4])
x2 = compose(M, ϵSE2, exp(M, ϵSE2, X2))

## ======================================================================================
## with Optim "point" as LieAlgebra coordinates and only one variable at a time as non-parametric
## ======================================================================================
using Manifolds
using Optim

M = SpecialEuclidean(2)

representation_size.(M.manifold.manifolds)


# algorithm=Optim.BFGS
algorithm=Optim.NelderMead 
alg = algorithm()
manifold = ManifoldsMani(SpecialEuclidean(2))
alg = algorithm(;manifold, algorithmkwargs...)

test_retract = false
function Optim.retract!(MM::ManifoldsMani, X)
    test_retract && (X[3] = rem2pi(X[3], RoundNearest))
    return X 
end

function Optim.project_tangent!(MM::ManifoldsMani, G, x)
    return G
end

options = Optim.Options(allow_f_increases=true, 
                        iterations = 200,
                        time_limit = 100,
                        show_trace = true,
                        show_every = 10,
                        )

# NOTE 
# Only for Lie Groups
# Input: Lie Algebra coordinates 
# f: exp to Group and calcuate residual
# Output: sum of squares residual 

function cost(X) 
    x = exp(M, ϵSE2, hat(M, ϵSE2, X))    
    return PriorPose2(priorx1, x)^2 + PriorPose2(priorx2, x)^2 + Pose2Pose2(measx1x2, x, x2)^2
end

initValues = @MVector [0.,0.,0.0]

@time result = Optim.optimize(cost, initValues, alg, options)

rv = Optim.minimizer(result)

## result on group:

rv_G = exp(M, ϵSE2, hat(M, ϵSE2, rv)) 

##
autodiff = :forward
# autodiff = :finite

initValues = [0.,0.,0.1]

tdtotalCost = Optim.TwiceDifferentiable((x)->cost(x), initValues; autodiff)

@time result = Optim.optimize(tdtotalCost, initValues, alg, options)

rv = Optim.minimizer(result)


##
H = Optim.hessian!(tdtotalCost, rv)
Σ = pinv(H)


## ======================================================================================
## with Optim flattened ProductRepr and only one variable at a time as non-parametric
## ======================================================================================
using Manifolds
using Optim

M = SpecialEuclidean(2)

representation_size.(M.manifold.manifolds)

struct ManifoldsMani <: Optim.Manifold
  mani::AbstractManifold
end


# flatten should be replaced by a view or @cast 
function flatten(M::SpecialEuclidean{2}, p)
    return mapreduce(vec, vcat, p.parts)
end

fp = flatten(M, p)

function unflatten(M::SpecialEuclidean{2}, fp::MVector)
    ProductRepr(MVector{2}(fp[1:2]), MMatrix{2,2}(fp[3:6]))
end

_p = unflatten(M, fp)


function Optim.retract!(MM::ManifoldsMani, fx)
    M = MM.mani
    x = unflatten(M, fx)
    project!(M, x, x)
    fx .= flatten(M, x)
    return fx 
end

function Optim.project_tangent!(MM::ManifoldsMani, fG, fx)
    M = MM.mani
    x = unflatten(M, fx)
    G = unflatten(M, fG)
    project!(M, G, x, G)
    fG .= flatten(M, G)
    return fG
end


# autodiff = :forward
algorithm=Optim.BFGS
autodiff = :finite
# algorithm=Optim.NelderMead # does not work with manifolds
algorithmkwargs=() # add manifold to overwrite computed one
options = Optim.Options(allow_f_increases=true, 
                        iterations = 100,
                        time_limit = 100,
                        show_trace = true,
                        show_every = 1,
                        )
##

function cost(fx) 
    x = unflatten(M, fx)
    return PriorPose2(priorx1, x)^2 + PriorPose2(priorx2, x)^2
end

manifold = ManifoldsMani(SpecialEuclidean(2))
alg = algorithm(;manifold, algorithmkwargs...)
# alg = algorithm(; algorithmkwargs...)


initValues = MVector(flatten(M, ϵSE2))
# tdtotalCost = Optim.TwiceDifferentiable((x)->cost(x), initValues; autodiff)
# result = Optim.optimize(tdtotalCost, initValues, alg, options)


result = Optim.optimize(cost, initValues, alg, options)

rv = Optim.minimizer(result)
unflatten(M, rv)


