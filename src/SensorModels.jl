# Linear array sonar constraints

# typealias FloatInt Union{Float64, Int64}

type LinearRangeBearingElevation <: Pairwise
  range::Normal
  bearing::Normal
  elev::Uniform
  LinearRangeBearingElevation() = new()
  LinearRangeBearingElevation( r::Tuple{Float64,Float64}, b::Tuple{Float64,Float64}; elev=Uniform(-0.25133,0.25133)) = new(Normal(r...),Normal(b...),elev)
end
function getSample( las::LinearRangeBearingElevation )
	return [rand(las.range); rand(las.bearing); rand(las.elev)]
end


# returns [Range Bearing Elevation] manifold difference between pose X ominus landmark L
function ominus(::Type{LinearRangeBearingElevation}, X::Vector{Float64}, L::Vector{Float64})
  # rangeBearing3(X, L)
  wTb = SE3(X[1:3], Euler(X[4:6]...))
  bTl = matrix(wTb)\[L[1:3];1.0]
  b = atan2(bTl[2],	bTl[1])
  el = -atan2(bTl[3], bTl[1])
  return [norm(bTl[1:3]); b; el]
end

# residual should equal zero when system is in balance
# measurement z is measurement vector with [range; bearing; elevation]
# variables are tuple (pose X [dim6], landmark L [dim3])
# function handle follows required parameter list
function residualLRBE!(residual::Vector{Float64}, z::Vector{Float64}, variables::Tuple)
  # TODO upgrade so the - sign here is used on a manifold too, ominus(z,  ominus(tt, variables...)  )
  residual[:] = z - ominus(LinearRangeBearingElevation, variables...)
  # @show residual
  nothing
end


# Convolution of conditional to project landmark location from position X (dim6)
# Y (dim3) are projected points
function project!(meas::LinearRangeBearingElevation, pose::Array{Float64,2}, landmark::Array{Float64,2}, idx::Int)
  landmark[1:3,idx] = numericRootGenericRandomized(residualLRBE!, getSample(meas), pose[1:6,idx][:], landmark[1:3,idx][:], 3) # ( bearrange3!,
	# x0[1:3,idx] = numericRoot(bearrange3!, getSample(meas), fixed[1:6,idx], x0[1:3,idx]+0.01*randn(3))
	nothing
end

function backprojectRandomized!(meas::LinearRangeBearingElevation, landmark::Array{Float64,2}, pose::Array{Float64,2}, idx::Int)
	pose[1:6,idx] = numericRootGenericRandomized(residualLRBE!, getSample(meas), landmark[1:3,idx][:], pose[1:6,idx][:], 3) # ( bearrange3!,
	nothing
end

# convenience function for project!
function project(meas::LinearRangeBearingElevation, X::Vector{Float64}, Y::Vector{Float64})
  xx, yy = X', Y'
  project!(meas, xx', yy', 1)
  return yy[1:3]
end

function evalPotential(meas::LinearRangeBearingElevation, Xi::Array{Graphs.ExVertex,1}, Xid::Int64)
  warn("evalPotential(meas::LinearRangeBearingElevation")
  fromX, ret, ff = zeros(0,0), zeros(0,0), +

  if Xi[1].index == Xid
    fromX = getVal( Xi[2] )
    ret = deepcopy(getVal( Xi[1] ))
    ff = backprojectRandomized!
  elseif Xi[2].index == Xid
    fromX = getVal( Xi[1] )
    ret = deepcopy(getVal( Xi[2] ))
    ff = project!
  end
  r,c = size(fromX)

	for i in 1:200
		ff(meas, fromX, ret, i)
	end

  return ret
end


## old code

# function backprojectRandomized!(meas::LinearRangeBearingElevation, landmark::Array{Float64,2}, pose::Array{Float64,2}, idx::Int)
# .....
# Zbr = sample(meas)
#
#   p = collect(1:6);
#   shuffle!(p);
#   p1 = p.==1; p2 = p.==2; p3 = p.==3
#   #@show x0, par
#   r = nlsolve(    (x, res) -> bearrangDidson!(res, Zbr,
#                   shuffleXAltD(x, x0[:,idx][:], 3, p), fixed[1:3,idx][:] ),
#                   [x0[p1,idx];x0[p2,idx];x0[p3,idx]]+0.01*randn(3)   )
#
#   X[1:6,idx] = shuffleXAltD(r.zero, x0[:,idx][:], 3, p );

# # randomized backproject using residual function
# function bearrange3!(residual::Vector{Float64}, Zrb::Vector{Float64}, X::Vector{Float64}, L::Vector{Float64})
#   wTb = SE3(X[1:3], Euler(X[4:6]...))
#   bTl = matrix(wTb)\[L[1:3];1.0]
#   b = atan2(bTl[2],	bTl[1])
# 	el = -atan2(bTl[3], bTl[1])
#   residual[1] = Zrb[1]-norm(bTl[1:3])
#   residual[2] = Zrb[2]-b
# 	residual[3] = Zrb[3]-el
#   nothing
# end

# returns vector with [range, bearing, elevation]
# takes pose X [dim6] and landmark L [dim3]
# function rangeBearing3( X::Vector{Float64}, L::Vector{Float64})
#   wTb = SE3(X[1:3], Euler(X[4:6]...))
#   bTl = matrix(wTb)\[L[1:3];1.0]
#   b = atan2(bTl[2],	bTl[1])
# 	el = -atan2(bTl[3], bTl[1])
#   return [norm(bTl[1:3]); b; el]
# end


# function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
# 		n = length(X)
#     Y = deepcopy(Alt)
# 		for i in 1:d
# 			Y[p[i]] = X[i]
# 		end
#     return Y
# end
#
# # use residual function to approximate the convolution of conditional belief with existing
# # belief estimate from fixed to x. Conditional belief is described by Pairwise measurement
# # zDim == length(sample(measurement))
# function numericRootGenericRandomized(residFnc::Function, measurement::IncrementalInference.Pairwise,
# 			fixed::Vector{Float64}, x0::Vector{Float64}, zDim::Int; perturb::Float64=0.01 )
#
# 	z = sample(meas)
# 	if zDim < length(x0)
# 		p = collect(1:6);
# 		shuffle!(p);
# 		# p1 = p.==1; p2 = p.==2; p3 = p.==3
# 		r = nlsolve(    (x, res) -> residFnc(res, Z,
#                     shuffleXAltD(x, x0, zDim, p), fixed ),
#                     x0[p[1:zDim]] + preturb*randn(zDim)   )
#     return shuffleXAltD(r.zero, x0, zDim, p );
# 	else
#     return (nlsolve(   (X, residual) -> residFnc(residual, measurement, parameters, X), x0 )).zero
# 	end
#   nothing
# end




# residual function must have the form
# bearrangDidson!(residual::Array{Float64,1}, Zbr::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
# function numericRoot(residFnc::Function, measurment, parameters, x0::Vector{Float64})
# 	return (nlsolve(   (X, residual) -> residFnc(residual, measurement, parameters, X), x0 )).zero
# end
# function solveLandmDidsonPt(measurement::Array{Float64,1}, param::Array{Float64,1}, x0::Array{Float64,1})
# 	numericRoot(bearrangDidson!, measurement, param, x0)
#   # return (nlsolve(   (X, residual) -> bearrangDidson!(residual, measurement, param, X), x0 )).zero
# end
