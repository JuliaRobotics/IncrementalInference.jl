
using Optim

# Normal(x; μ, σ) = 1/(σ/sqrt(2π)) exp(-0.5((μ-x)/σ)^2 )
#
# two Normals (z1, μ1=-1, z2, μ2=+1) on one variable (x)
# z1 ⟂ z2

# [Θ | z] ∝ [z1 | x] × [z2 | x]
# [Θ | z] ∝ N(x;-1, 1) × N(x; +1, 1)
# N(x; μ, σ) ∝ N(x; 0, 1/sqrt(2))
#
# 1/σ^2 = 1/σ1^2 + 1/σ2^2
# μ/σ^2 = μ1/σ1^2 + μ2/σ2^2
#

# Malahanobis distance
# res' * iΣ * res
# [δx]' * [1/sigma^2] * [δx]

# - log
f(x) = (x[1]-1)^2 + (x[1]+1)^2
# f(x) = [x-1]'*1/1^2*[x-1] + [x+1]'*1/1^2*[x+1]

f(x) = 1/1*(x[1])^2      # 0.5   1/2 × 2 = σ
f(x) = 1/0.5*(x[1])^2    # 0.25  1/4 × 2 = σ
f(x) = 1/0.25*(x[1])^2   # 0.125 1/8 × 2 = σ

init = [0.0]

func = TwiceDifferentiable(f, init)#; autodiff=:forward);

opt = optimize(func, init)

parameters = Optim.minimizer(opt)

numerical_hessian = hessian!(func, parameters)

cov_matrix = pinv(numerical_hessian) # inv gives the same answer


## Sample one Guassian...

using Optim, NLSolversBase, Random
using LinearAlgebra: diag
Random.seed!(0);                            # Fix random seed generator for reproducibility

n = 500                             # Number of observations
nvar = 1                            # Number of variables
β = ones(nvar) * 3.0                # True coefficients
x = [ones(n) randn(n, nvar - 1)]    # X matrix of explanatory variables plus constant
ε = randn(n) * 0.5                  # Error variance
y = x * β + ε;                      # Generate Data

function Log_Likelihood(X, Y, β, log_σ)
    σ = exp(log_σ)
    llike = -n/2*log(2π) - n/2* log(σ^2) - (sum((Y - X * β).^2) / (2σ^2))
    llike = -llike
end

function Log_Likelihood(X, μ, σ)
  log(  1/(σ/sqrt(2π)) * exp(-0.5(μ-X)/σ)^2  )
end



func = TwiceDifferentiable(vars -> Log_Likelihood(x, y, vars[1:nvar], vars[nvar + 1]),
                           ones(nvar+1); autodiff=:forward);

opt = optimize(func, ones(nvar+1))

parameters = Optim.minimizer(opt)

parameters[nvar+1] = exp(parameters[nvar+1])

numerical_hessian = hessian!(func,parameters)

var_cov_matrix = inv(numerical_hessian)

β = parameters[1:nvar]

temp = diag(var_cov_matrix)
temp1 = temp[1:nvar]

t_stats = β./sqrt.(temp1)




## two gaussians in different dimensions
# http://mathworld.wolfram.com/NormalProductDistribution.html
# Bessel K0
# using SpecialFunctions
# besselk(0, 1)
