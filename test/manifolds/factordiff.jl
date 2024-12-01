
# dev test script for factor gradients using ForwardDiff and ManifoldDiff

using ManifoldDiff
using Manifolds
using IncrementalInference
using StaticArrays
using Zygote
using Test

##


@testset "test factorJacobian" begin
##

# manual EuclidDistance 
z = 10.0
f_eucliddist(x1,x2; z) = z - norm(x2 - x1)
x0, x1 = [0.0; 0.0], [10.0, 0.0]
J_ = Zygote.jacobian(()->f_eucliddist(x0, x1; z), Zygote.Params([x0, x1]))

Jx0 = J_[x0]
Jx1 = J_[x1]


##

fg = LocalDFG(;
  solverParams = SolverParams(;
    graphinit=false
  )
)

addVariable!.(fg, [:x0; :x1], Position2)
f = addFactor!(fg, [:x0; :x1], EuclidDistance(Normal(z,1.0)))

p1 = [SA[x1[1];x1[2]] for _ in 1:1]

setVal!(fg, :x1, p1, solveKey=:parametric)

J = IIF.factorJacobian(fg, :x0x1f1)

@test isapprox( Jx0, J[1:1,1:2]; atol=1e-8)
@test isapprox( Jx1, J[1:1,3:4]; atol=1e-8)


##
end
##


# ## 
# @testset "using RoME; FiniteDiff.jacobian of SpecialEuclidean(2) factor" begin
# ##

# fg = LocalDFG(;
#   solverParams = SolverParams(;
#     graphinit=false
#   )
# )

# addVariable!.(fg, [:x0; :x1], Pose2)
# f = addFactor!(fg, [:x0; :x1], Pose2Pose2(MvNormal([10;0;pi/2],[1 0 0; 0 1 0; 0 0 1.0])))

# p1 = [ArrayPartition([10; 0.0], [0 1; -1 0.0]) for _ in 1:1]

# setVal!(fg, :x1, p1, solveKey=:parametric)

# J = IIF.factorJacobian(fg, :x0x1f1)

# @test isapprox( Jx0, J[1:1,1:2]; atol=1e-8)
# @test_broken isapprox( Jx1, J[1:1,3:4]; atol=1e-8)


# ##
# end
# ##


## 
@testset "ManifoldDiff.jacobian of SpecialEuclidean(2) factor" begin
##


M = SpecialEuclidean(2; vectors=HybridTangentRepresentation())
z = ArrayPartition(SA[10.0; 0.0], SMatrix{2,2}(0.0, -1.0, 1.0, 0.0))

p1 = ArrayPartition(SA[0.0; 0.0], SMatrix{2,2}(1, 0, 0, 1.))
e0 = identity_element(M, p1)
p2 = exp(M, e0, hat(M, e0, [10,0,pi/2]))


function resid_SE2(X, p, q)
  q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
  return vee(M, q, log(M, q, q̂))
end


# finitediff setup
# finitediff setup
r_backend = ManifoldDiff.TangentDiffBackend(
  if v"0.4" <=  pkgversion(ManifoldDiff)
    ManifoldDiff.AutoFiniteDifferences(central_fdm(5, 1))
  else
    ManifoldDiff.FiniteDifferencesBackend()
  end
)
Me = Euclidean(3)

function _factorJac!(J, z, p1, p2)
  g1(p_) = resid_SE2(z, p_, p2)
  g2(p_) = resid_SE2(z, p1, p_)
  
  J[1:3,1:3] = ManifoldDiff.jacobian(M, Me, g1, p1, r_backend)
  J[1:3,4:6] = ManifoldDiff.jacobian(M, Me, g2, p2, r_backend)
  
  J
end
# f_SE2_z(z_) = resid_SE2(z_, e0, p2)
# f_SE2_x0(p_) = resid_SE2(z, e0, p_)
# f_SE2_x0(p_) = resid_SE2(z, e0, p_)

J = zeros(3,6)

J_ = _factorJac!(J, z, p1, p2)
# @profview _factorJac!(J, z, p1, p2)

if false
  # finitediff setup
  z_backend = ManifoldDiff.TangentDiffBackend(
    if v"0.4" <=  pkgversion(ManifoldDiff)
      ManifoldDiff.AutoFiniteDifferences(central_fdm(5, 1))
    else
      ManifoldDiff.FiniteDifferencesBackend()
    end
  )

  g = ManifoldDiff.jacobian(M, Euclidean(3), f_SE2_x0, p1, z_backend)
else
  @info "ManifoldDiff.ZygoteDiffBackend usage still under development (23Q3)"
end



##
end

# cost(p_) = distance(M, e0, p_) # ManifoldDiff.gradient(M, cost, p, r_backend)
# cost(p)
# ManifoldDiff.gradient(M, cost, p, r_backend)










# ##
# @testset "CODE STILL UNDER DEV:::Zygote on SpecialEuclidean(2)" begin
# ##

# # manual PosePose 
# M = SpecialEuclidean(2)
# z = ArrayPartition(SA[10.0; 0.0], SMatrix{2,2}(0.0, -1.0, 1.0, 0.0))
# e0 = identity_element(M, z)

# # modified from IIF/test/testSpecialEuclidean2Mani.jl
# function f_SE2(X, p, q)
#   q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
#   Xc = zeros(3)
#   vee!(M, Xc, q, log(M, q, q̂))
#   return Xc
#   # return (Xc'*Xc)[1]
# end

# Xc_0 = [0.0; 0.0; 0.0] # deepcopy(e0)
# Xc_1 = [10.0; 0.0; pi/2] # deepcopy(z)

# J_ = Zygote.jacobian(
#   ()->begin
#     f_SE2(
#       log(M, e0, z), 
#       exp(M, e0, hat(M, e0, Xc_0)), 
#       exp(M, e0, hat(M, e0, Xc_1))
#     )
#   end, 
#   Zygote.Params([Xc_0, Xc_1])
# )

# Jx0 = J_[x0]
# Jx1 = J_[x1]


# ##
# end