# from #605

# using Revise

using IncrementalInference
using Test
using TensorCast

##

@testset "test mixture prior" begin
##

# init graph
fg = initfg()
N = 100
getSolverParams(fg).N = N;

# add first variable
addVariable!(fg, :x0, ContinuousScalar)

# also test AliasingScalingSampler
v = rand(50)
v[20:29] .+= 5*rand(10)
v ./= sum(v)
bss = AliasingScalarSampler(collect(1:50), v)


# add bi-modal mixture prior
Prior0 = Mixture(Prior,(a=Normal(-5.0,1.0), b=Uniform(0.0,1.0)), [0.5,0.5])
Prior0 = Mixture(Prior,(Normal(-5.0,1.0), Normal(0.0,1.0)), [0.5,0.5])
Prior0 = Mixture(Prior,[Normal(-5.0,1.0), Normal(0.0,1.0)], [0.5,0.5])
Prior0 = Mixture(Prior,(Normal(-5.0,1.0), Normal(0.0,1.0)), [0.5;0.5])
Prior0 = Mixture(Prior,[Normal(-5.0,1.0), Normal(0.0,1.0)], [0.5;0.5])
Prior0 = Mixture(Prior,(Normal(-5.0,1.0), bss), Categorical([0.5;0.5]))
f1 = addFactor!(fg, [:x0], Prior0)

# also test serialization of AliasingScalarSampler
saveDFG("/tmp/test_fg_bss", fg)


# check numerics -- replaced by CalcFactor{<:Mixture}, 
smpls_ = approxConv(fg, :x0f1, :x0)
# smpls, = getSample(Prior0, N) # ,lb=

# lazy
@cast smpls[i,j] := smpls_[j][i]


# should be a balance of particles
# @test sum(lb .== 1) - sum(lb .== 2) |> abs < 0.3*N
@test sum(smpls .< -2.5) - sum(-2.5 .< smpls) |> abs < 0.35*N

# solve
solveTree!(fg);

marginalPts_ = getBelief(fg, :x0) |> getPoints

# lazy
@cast marginalPts[i,j] := marginalPts_[j][i]

# check solver solution consistent too
@test sum(marginalPts .< -2.5) - sum(-2.5 .< marginalPts) |> abs < 0.35*N

##
end


@testset "Serialization of Mixture(Prior,..) including a AliasingScalarSampler" begin
##

fg_ = loadDFG("/tmp/test_fg_bss")

N = getSolverParams(fg_).N


solveTree!(fg_);

marginalPts_ = getBelief(fg_, :x0) |> getPoints

# lazy
@cast marginalPts[i,j] := marginalPts_[j][i]

# check solver solution consistent too
@test sum(marginalPts .< -2.5) - sum(-2.5 .< marginalPts) |> abs < 0.35*N


# cleanup
Base.rm("/tmp/test_fg_bss.tar.gz")

##
end


# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,20cm)

# # plot the results
# plotKDE(fg, :x0)


#
