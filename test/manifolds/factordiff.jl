
# dev test script for factor gradients using ForwardDiff and ManifoldDiff

# using ManifoldDiff
# using Manifolds
using IncrementalInference
using StaticArrays
# using DataStructures: OrderedDict

##

# M = SpecialEuclidean(2)
# p = ArrayPartition(SA[0.0; 0.0], SMatrix{2,2}(1, 0, 0, 1.))

# # finitediff setup
# r_backend = ManifoldDiff.TangentDiffBackend(
#   ManifoldDiff.FiniteDifferencesBackend()
# )

# e0 = identity_element(M, p)

# cost(p_) = distance(M, e0, p_)

# cost(p)

# g = ManifoldDiff.jacobian(M, TranslationGroup(3), cost, p, r_backend)


##

fg = LocalDFG(;
  solverParams = SolverParams(;
    graphinit=false
  )
)

addVariable!.(fg, [:x0; :x1], Position2)
f = addFactor!(fg, [:x0; :x1], EuclidDistance(Normal(10,1.0)))

p1 = [SA[0.0;0.0] for _ in 1:1]

setVal!(fg, :x1, p1, solveKey=:parametric)


J = factorJacobian(fg, :x0x1f1)

##
