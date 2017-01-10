# try threads
using Base.Test

@show Threads.nthreads()

function ff!(ret::Vector{Float64}, x::Vector{Float64}, i::Int)
  ret[i] = x[i]
  nothing
end

N = 12
a = randn(N)
r = zeros(N)
for i in 1:N  ff!(r,a,i)  end

r = zeros(N)
@time for i in 1:N  ff!(r,a,i)  end
r = zeros(N)
@time for i in 1:N  ff!(r,a,i)  end
@test norm(r-a) < 1e-10

r = zeros(N)
Threads.@threads for i in 1:N  ff!(r,a,i)  end

r = zeros(N)
@time Threads.@threads for i in 1:N  ff!(r,a,i)  end
@test norm(r-a) < 1e-10


r = zeros(N)
@time for i in 1:N
  r[i] = a[i]^2
end


r = zeros(N)
@time Threads.@threads for i in 1:N
  r[i] = a[i]^2
end




using IncrementalInference, TransformUtils

# see speedup inpact here
# 3 dimensional line, z = [a b][x y]' + c
function rotationresidual!(res::Vector{Float64}, z::Vector{Float64}, var::Tuple)
  q1 = convert(Quaternion, Euler(z...))
  q2 = convert(Quaternion, so3(var[2]))
  qq = q1*q_conj(q2)
  res[1:3] = vee(convert(so3, qq))
  nothing
end

N=12
eul = 0.25*randn(3,N)
y = zeros(3,N)

@time for i in 1:N
  gg = (x, res) -> rotationresidual!(res, eul[:,i], (zeros(0),x))
  y[:,i] = numericRootGenericRandomizedFnc(
          gg,
          3, 3, 0.1*randn(3)    )
  # test the result
  # @test TransformUtils.compare(convert(Quaternion, Euler(eul...)),
  #                             convert(Quaternion, so3(y)), tol=1e-8)
end












#
