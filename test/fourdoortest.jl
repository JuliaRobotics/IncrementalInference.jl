using IncrementalInference
using Test
##

@testset "test fourdoor early example" begin

## example parameters
# Number of kernels representing each marginal belief
N=100

# prior knowledge of four possible door locations
cv = 3.0
doorPrior = Mixture(Prior,
                    [Normal(-100,cv);Normal(0,cv);Normal(100,cv);Normal(300,cv)],
                    [1/4;1/4;1/4;1/4] )

## Build the factor graph object
fg = initfg()
getSolverParams(fg).useMsgLikelihoods = true

# first pose location
v1 = addVariable!(fg,:x1,ContinuousScalar,N=N)
# see a door for the first time
addFactor!(fg,[:x1], doorPrior)

# first solution with only one variable and factor (may take a few moments on first JIT compiling)
solveTree!(fg)


## drive to second pose location

addVariable!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[:x1;:x2],LinearRelative(Normal(50.0,2.0)))

# drive to third pose location
v3=addVariable!(fg,:x3,ContinuousScalar, N=N)
addFactor!(fg,[:x2;:x3], LinearRelative( Normal(50.0,4.0)))

# see a door for the second time
addFactor!(fg,[:x3], doorPrior)

# second solution should be much quicker
solveTree!(fg)

# drive to forth and final pose location
addVariable!(fg,:x4,ContinuousScalar, N=N)
addFactor!(fg,[:x3;:x4], LinearRelative( Normal(200.0,4.0)))


## make a third door sighting
addFactor!(fg,[:x4], doorPrior)

# solve over all data
tree = solveTree!(fg)

##


end



# # HMM computed ground truth, extended for 7 poses with landmark
# global gt = Dict{Symbol, Array{Float64,2}}()
# gt[:x0]=reshape(Float64[0.0;1.97304 ],2,1) # -0.0342366
# gt[:x2]=reshape(Float64[50.0; 2.83153 ],2,1) # 49.8797
# gt[:x3]=reshape(Float64[100.0; 1.65557 ],2,1) # 99.8351
# gt[:x4]=reshape(Float64[150.0; 1.64945 ],2,1) # 148.637
# gt[:x5]=reshape(Float64[200.0; 1.77992 ],2,1) # 198.62
# gt[:x6]=reshape(Float64[240.0; 2.20466 ],2,1) # 238.492
# gt[:x7]=reshape(Float64[300.0; 2.14353 ],2,1) # 298.467
# gt[:l1]=reshape(Float64[165.0; 1.17284 ],2,1) # 164.102
#
