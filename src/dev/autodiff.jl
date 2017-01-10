# test autodiff

using ForwardDiff
using TransformUtils

# pose3pose3
x = zeros(2)

cfg = ForwardDiff.JacobianConfig{2}(x)

function ff(x::Vector)
  return [sin(x[1]);cos(x[2])]
end

ForwardDiff.jacobian(ff, x, cfg)


function ff2(res::Vector{Float64}, idx::Int, X::Array{Float64,2})
  res[1] = sin(X[1,idx]);
  res[2] = cos(X[2,idx]);
  return
end

type FF
  ret::Vector{Float64}
  idx::Int
  X::Array{Float64,2}
  usrfnc!::Function
end

function (fd::FF)(x::Vector)
  fd.X[:,fd.idx] = x[:]
  fd.usrfnc!(fd.ret, fd.idx, fd.X)
  return fd.ret
end

ffvar = FF(zeros(2), 1, zeros(2,3), ff2)

ffvar(x)
ForwardDiff.jacobian(ffvar, x, cfg) # this does not work


using NLsolve

function f!(x, fvec)
    fvec[1] = (x[1]+3)*(x[2]^3-7)+18
    fvec[2] = sin(x[2]*exp(x[1])-1)
end

function g!(x, fjac)
    fjac[1, 1] = x[2]^3-7
    fjac[1, 2] = 3*x[2]^2*(x[1]+3)
    u = exp(x[1])*cos(x[2]*exp(x[1])-1)
    fjac[2, 1] = x[2]*u
    fjac[2, 2] = u
end

@time nlsolve(f!, [ 0.1; 1.2], autodiff=true)
