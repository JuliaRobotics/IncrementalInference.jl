using IncrementalInference
using IncrementalInference.Optim
using Test

##
@testset "test parametric mixture" begin

## Test simple mixture prior

fg = LocalDFG(;solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x0, ContinuousScalar)

p = Mixture(Prior(I), [Normal(0.8, 0.4), Normal(1.0, 0.1)], Categorical([0.1; 0.9]))
f = addFactor!(fg, [:x0], p)

addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0,:x1], LinearRelative(Normal(1.0,0.1)))

options = Optim.Options(time_limit = 100,
                    iterations = 10000,
                    show_trace = true,
                    show_every = 10,
                    allow_f_increases=true,
                    g_tol = 1e-6,
                    )

algorithm = Optim.NelderMead
vardict, result, varIds, Σ = IIF.solveGraphParametric(fg; options, algorithm)

@test isapprox(vardict[:x0].val[1], 1, atol = 0.05)
@test isapprox(vardict[:x0].cov[1], 0.01, atol = 0.001)


##
fg = LocalDFG(;solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal(0.0,0.1)))

mlr = Mixture(LinearRelative(I), [Normal(-1.0, 0.2), Normal(1.0, 0.1)], Categorical([0.5; 0.5]))

addFactor!(fg, [:x0,:x1], mlr)

addVariable!(fg, :l1, ContinuousScalar)
p = Mixture(Prior(I), [Normal(-1.5, 0.1), Normal(0.9, 0.2)], Categorical([0.5; 0.5]))
# p = Prior(Normal(0.9, 0.1))
addFactor!(fg, [:l1], p)

addFactor!(fg, [:x1,:l1], LinearRelative(Normal(0.0,0.1)))

vardict, result, varIds, Σ = IIF.solveGraphParametric(fg)
@test isapprox(vardict[:x0].val[1], 0, atol = 0.1)
@test isapprox(vardict[:x1].val[1], 1, atol = 0.1)
@test isapprox(vardict[:l1].val[1], 1, atol = 0.1)


## ContinuousEuclid(2) prior

fg = LocalDFG(;solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x0, ContinuousEuclid(2))

p = Mixture(Prior(MvNormal(2,1.0)), [MvNormal([0.8, 0.5], [0.4, 0.4]), MvNormal([1.0, 0.5], [0.1, 0.1])], Categorical([0.1; 0.9]))
f = addFactor!(fg, [:x0], p)


vardict, result, varIds, Σ = IIF.solveGraphParametric(fg)
vardict
@test isapprox(vardict[:x0].val[1], 1.0, atol = 0.01)
@test isapprox(vardict[:x0].val[2], 0.5, atol = 0.01)
@test isapprox(vardict[:x0].cov[1,1], 0.01, atol = 0.001)
@test isapprox(vardict[:x0].cov[2,2], 0.01, atol = 0.001)

# x = collect(0.5:0.01:1.5)
# y = collect(0.0:0.01:1.0)
# cfm = IIF.CalcFactorMahalanobis(f)
# xy = Vector{Float64}[]
# x = 0.5:0.01:1.5
# y = 0:0.01:1
# for x=x,y=y
#     push!(xy, [x,y]) 
# end
# r = cfm.(xy)
# surface(x,y,r)

##
end











if false

fg = LocalDFG(;solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal(0.0,0.1)))

mlr = Mixture(LinearRelative(I), [Normal(-1.0, 0.2), Normal(1.0, 0.1)], Categorical([0.5; 0.5]))
# mlr = Mixture(LinearRelative(I), [Normal(-0.2, 0.05), Normal(0.2, 0.1)], Categorical([0.5; 0.5]))
# mlr = Mixture(LinearRelative(I), [Normal(0.0, 0.9), Normal(0.0, 0.1)], Categorical([0.5; 0.5]))

# mlr = Mixture(LinearRelative(I), [Normal(-0.2, 0.1), Normal(0.2, 0.1)], Categorical([0.9; 0.1]))
# mlr = Mixture(LinearRelative(I), [Normal(-0.2, 0.15), Normal(0.2, 0.1)], Categorical([0.5; 0.5]))
# mlr = Mixture(LinearRelative(I), [Normal(-1.0, 0.2), Normal(1.0, 0.1)], Categorical([0.5; 0.5]))

addFactor!(fg, [:x0,:x1], mlr)

f = fg[:x0x1f1]
cfm = IIF.CalcFactorMahalanobis(f)
x1 = collect(-2:0.01:2)
# x1 = reverse(collect(-2:0.01:2))
# x1 = collect(-0.5:0.01:0.5)
r1 = cfm.(Ref([0.0]), x1)
# plot!(x1,r)
plot(x1,r1)

##
addVariable!(fg, :l1, ContinuousScalar)
p = Mixture(Prior(I), [Normal(-1.5, 0.1), Normal(0.9, 0.2)], Categorical([0.5; 0.5]))
# p = Prior(Normal(0.9, 0.1))
addFactor!(fg, [:l1], p)

addFactor!(fg, [:x1,:l1], LinearRelative(Normal(0.0,0.1)))

f = fg[:l1f1]
cfm = IIF.CalcFactorMahalanobis(f)
x1 = collect(-2:0.01:2)
r2 = cfm.(x1)
plot!(x1,r2)
# plot!(x1,r1+r2)

##
initAll!(fg)
solveTree!(fg)
plotKDE(fg, :x1)
plotKDE(fg, :l1)
##
using IncrementalInference.Optim

options = Optim.Options(time_limit = 100,
                    iterations = 1000,
                    show_trace = true,
                    show_every = 1,
                    allow_f_increases=true,
                    g_tol = 1e-6,
                    )

algorithm = Optim.NelderMead
algorithm = Optim.BFGS
algorithmkwargs=(linesearch = Optim.LineSearches.HagerZhang(),)
# algorithmkwargs=(linesearch = Optim.LineSearches.Static(),)
vardict, result, varIds, Σ = IIF.solveGraphParametric(fg; algorithm, options, algorithmkwargs)
vardict

##

IIF.updateParametricSolution!(fg, vardict)

##

end

## 
if false
using RoME

fg = LocalDFG(;solverParams=SolverParams(algorithms=[:default, :parametric]))

pr_noise = [0.01, 0.01, 0.001]
od_noise_1 = [0.5; 0.5; 0.2]
od_noise_2 = [0.05; 0.05; 0.02]

#x0 prior
addVariable!(fg, :x0, Pose2)
prpo = MvNormal([0,0,-pi], pr_noise)
addFactor!(fg, [:x0], PriorPose2(prpo))
# addFactor!(fg, [:x0], PriorPose2(MvNormal(rand(prpo), pr_noise)))

#x0 to x1
addVariable!(fg, :x1, Pose2)
pp1 = MvNormal([1.0,0,0], od_noise_1)
pp2 = MvNormal([1.0,0,0], od_noise_2)
mpp = Mixture(Pose2Pose2(I), [pp1, pp2], Categorical([0.5; 0.5]))

f = addFactor!(fg, [:x0,:x1], mpp)

cfm = IIF.CalcFactorMahalanobis(f)


##
xyθ = Vector{Float64}[]
for x = 0.5:0.01:1.5, y = 0, θ = 0
    push!(xyθ, [x,y,θ]) 
end

x = 0.5:0.01:1.5

r = cfm.(Ref([0.0; 0.0; 0.0]), xyθ)
plot(x,r)

##
#x1 to x2
# addVariable!(fg, :x2, Pose2)
# pp1 = MvNormal([1.0,0,0], od_noise_1)
# pp2 = MvNormal([1.0,0,0], od_noise_2)
# mpp = Mixture(Pose2Pose2(I), [pp1, pp2], Categorical([0.5; 0.5]))
# addFactor!(fg, [:x1,:x2], mpp)

# vardict, result, varIds, Σ = solveFactorGraphParametric(fg)
# vardict

##

using IncrementalInference.Optim

options = Optim.Options(time_limit = 100,
                    iterations = 10000,
                    show_trace = true,
                    show_every = 1,
                    allow_f_increases=true,
                    g_tol = 1e-6,
                    )

algorithm = Optim.NelderMead
vardict, result, varIds, Σ = IIF.solveGraphParametric(fg; algorithm, autodiff=:finite)
vardict


end
