# development testing for multithreaded convolution

using IncrementalInference
using Test
using TensorCast

##

@testset "Basic ContinuousScalar example to ensure multithreaded convolutions work..." begin

##

@show Threads.nthreads()

N = 100

# Start with an empty factor graph
fg = initfg()

# add the first node
addVariable!(fg, :x0, ContinuousScalar, N=N)

# this is unary (prior) factor and does not immediately trigger autoinit of :x0.
addFactor!(fg, [:x0], Prior(Normal(0,1)))


addVariable!(fg, :x1, ContinuousScalar, N=N)
# P(Z | :x1 - :x0 ) where Z ~ Normal(10,1)
@warn "threadmodel=MultiThreaded is obsolete.  Look at IIF.CalcFactor alternatives instead"
addFactor!(fg, [:x0, :x1], LinearRelative(Normal(10.0,1)) ) #, threadmodel=MultiThreaded)

@error "Factor threadmodel=MultiThreaded equivalence restoration TBD"
@test begin
    pts_ = approxConv(fg, :x0x1f1, :x1, N=N)
    @cast pts[i,j] := pts_[j][i]

    @test 0.95*N <= sum(abs.(pts .- 10.0) .< 5.0)
    true
end
##

end


#
