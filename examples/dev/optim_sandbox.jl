include("factors_sandbox.jl")

using Optim

# TODO 2 variables x1 and x2 of type SE2 one factor and one prior on x1

X = hat(M, ϵSE2, [1, 0, pi/4])
priorx1 = compose(M, ϵSE2, exp(M, ϵSE2, X))

X = hat(M, ϵSE2, [0, 0, pi/4])
priorx2 = compose(M, ϵSE2, exp(M, ϵSE2, X))

measX = hat(M, ϵSE2, [1, 0, pi/4])
measx1x2 = compose(M, ϵSE2, exp(M, ϵSE2, measX))


## ======================================================================================
## with Optim "point" as LieAlgebra coordinates and only one variable at a time as non-parametric
## ======================================================================================
using Manifolds
using Optim

M = SpecialEuclidean(2)

representation_size.(M.manifold.manifolds)

# autodiff = :forward
algorithm=Optim.BFGS
autodiff = :finite
# algorithm=Optim.NelderMead 
alg = algorithm()
options = Optim.Options(allow_f_increases=true, 
                        iterations = 100,
                        time_limit = 100,
                        show_trace = true,
                        show_every = 1,
                        )

# NOTE 
# Die gaan net op groepe werk
# tangent coordinates kom in
# gaan af na die manifold bereken en dan weer terug

function cost(X) 
    # @info X
    x = exp(M, ϵSE2, hat(M, ϵSE2, X))    
    return PriorPose2(priorx1, x)^2 + PriorPose2(priorx2, x)^2
end

initValues = @MVector [0.,0.,0.]

result = Optim.optimize(cost, initValues, alg, options)

rv = Optim.minimizer(result)


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

