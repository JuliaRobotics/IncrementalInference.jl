
function numericRoot(residFnc::Function, measurement, parameters, x0::Vector{Float64})
	return (nlsolve(   (X, res) -> residFnc(res, measurement, parameters, X), x0 )).zero
end


function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
		# n = length(X)
    Y = deepcopy(Alt)
    for i in 1:d
			Y[p[i]] = X[i]
		end
    return Y
end

# residual function must have the form residFnc!(res, x)
function numericRootGenericRandomizedFncOld(
      residFnc!::Function,
      zDim::Int,
      xDim::Int,
      x0::Vector{Float64};
      perturb::Float64=1e-5  )

	# xDim = length(x0)
	if zDim < xDim
		p = collect(1:xDim);
		shuffle!(p);
		# p1 = p.==1; p2 = p.==2; p3 = p.==3
    # TODO -- refactor to use functors instead
		r = nlsolve(    (x, res) -> residFnc!(shuffleXAltD(x, x0, zDim, p), res ),
                    x0[p[1:zDim]] + perturb*randn(zDim)
               )

    return shuffleXAltD(r.zero, x0, zDim, p );
	else
    # return (nlsolve(  (x, res) -> residFnc!(res, x),
    #                   x0      )).zero
    # @show x0
    return ( nlsolve(  residFnc!, x0 + perturb*randn(zDim) ) ).zero
	end
end

#-------------------------------------------------------------------------------
# first in place attempt for generic wrap and root finding below

# TODO -- <: FunctorPairwise seems to be a problem here, abstraction feels a bit wrong
type GenericWrapParam{T} <: FunctorPairwise
  usrfnc!::T
  params::Vector{Array{Float64,2}}
  varidx::Int
  particleidx::Int
  measurement::Array{Float64,2}
  samplerfnc::Function
  GenericWrapParam() = new()
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}) = new(fnc, t, 1,1, zeros(0,1), +)
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int) = new(fnc, t, i, j, zeros(0,1), +)
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Array{Float64,2}, smpl::Function) = new(fnc, t, i, j, meas, smpl)
end
function (p::GenericWrapParam)(x, res)
  # approximates by not considering cross indices among parameters
  p.params[p.varidx][:, p.particleidx] = x
  p.usrfnc!(res, p.particleidx, p.measurement, p.params...)
end

type FastRootGenericWrapParam{T} <: Function
  p::Vector{Int}
  perturb::Vector{Float64}
  X::Array{Float64,2}
  Y::Vector{Float64}
  xDim::Int
  zDim::Int
  gwp::GenericWrapParam{T}
  FastRootGenericWrapParam{T}(xArr::Array{Float64,2}, zDim::Int, residfnc::T) =
      new(collect(1:size(xArr,1)), zeros(zDim), xArr, zeros(size(xArr,1)), size(xArr,1), zDim, residfnc)
end

# Shuffle incoming X into random permutation in fr.Y
# shuffled fr.Y will be placed back into fr.X[:,fr.gwp.particleidx] upon fr.gwp.usrfnc(x, res)
function shuffleXAltD!(fr::FastRootGenericWrapParam, X::Vector{Float64})
  copy!(fr.Y, fr.X[:,fr.gwp.particleidx])
  for i in 1:fr.zDim
    fr.Y[fr.p[i]] = X[i]
  end
  nothing
end
function (fr::FastRootGenericWrapParam)( x::Vector{Float64}, res::Vector{Float64} )
  shuffleXAltD!(fr, x)
  fr.usrfnc( fr.Y, res ) #function (p::GenericWrapParam)(x, res)
end

# Solve free variable x by root finding residual function fgr.usrfnc(x, res)
# randomly shuffle x dimensions if underconstrained by measurement z dimensions
# small random perturbation used to prevent trivial solver cases, div by 0 etc.
# result stored in fgr.Y
# fr.X must be set to memory ref the param[varidx] being solved, at creation of fr
function numericRootGenericRandomizedFnc!{T}(
      fr::FastRootGenericWrapParam{T};
      perturb::Float64=1e-5,
      testshuffle::Bool=false )
  #
  fr.perturb[1:fr.zDim] = perturb*randn(fr.zDim)
	if fr.zDim < fr.xDim || testshuffle
    shuffle!(fr.p)
    r = nlsolve(  fr,
                  fr.X[fr.p[1:fr.zDim], fr.gwp.particleidx] + fr.perturb # this is x0
               )
    shuffleXAltD!( fr, r.zero )
	else
    fr.Y = ( nlsolve(  fr.gwp, fr.X[:,fr.gwp.particleidx] + fr.perturb ) ).zero
	end
  fr.X[:,fr.gwp.particleidx] = fr.Y
  nothing
end


#-------------------------------------------------------------------------------


# this will likely expand with more internal bells and whistles
# to perform in place memory operations for array values in

# see FastRootGenericWrapParam{T}

# TODO -- part of speed and flexibility refactoring exercise
type FastGenericRoot{T} <: Function
  p::Vector{Int}
  perturb::Vector{Float64}
  X::Vector{Float64}
  Y::Vector{Float64}
  xDim::Int
  zDim::Int
  usrfnc::T
  FastGenericRoot{T}(xDim::Int, zDim::Int, residfnc::T) =
      new(collect(1:xDim), zeros(zDim), zeros(xDim), zeros(xDim), xDim, zDim, residfnc)
end
function shuffleXAltD!(fr::FastGenericRoot, X::Vector{Float64})
  copy!(fr.Y, fr.X)
  for i in 1:fr.zDim
    fr.Y[fr.p[i]] = X[i]
  end
  nothing
end
function (fr::FastGenericRoot)( x::Vector{Float64}, res::Vector{Float64} )
  shuffleXAltD!(fr, x)  #(fr.Y, x, fr.x0, fr.zDim, fr.p)
  fr.usrfnc( fr.Y, res )
end
# Solve free variable x by root finding residual function fgr.usrfnc(x, res)
# randomly shuffle x dimensions if underconstrained by measurement z dimensions
# small random perturbation used to prevent trivial solver cases, div by 0 etc.
# result stored in fgr.Y
function numericRootGenericRandomizedFnc!{T}(
      fgr::FastGenericRoot{T};
      perturb::Float64=1e-5,
      testshuffle::Bool=false )
  #
  fgr.perturb[1:fgr.zDim] = perturb*randn(fgr.zDim)
	if fgr.zDim < fgr.xDim || testshuffle
    shuffle!(fgr.p)
    r = nlsolve(  fgr,
                  fgr.X[fgr.p[1:fgr.zDim]] + fgr.perturb # this is x0
               )
    # copy!(fgr.X, x0) # should need this line?
    shuffleXAltD!( fgr, r.zero )
    # copy!(fgr.X, r.zero)
	else
    fgr.Y = ( nlsolve(  fgr.usrfnc, fgr.X + fgr.perturb ) ).zero
    # copy!(fgr.X, y)
	end
  nothing
end


function numericRootGenericRandomizedFnc(
      residFnc!::Function,
      zDim::Int,
      xDim::Int,
      x0::Vector{Float64};
      perturb::Float64=1e-5,
      testshuffle::Bool=false   )
  #
  # TODO -- this only start of refactoring for inplace, more to come
  # xDim = length(x0)
  fgr = FastGenericRoot{typeof(residFnc!)}(xDim, zDim, residFnc!)
  shuffle!(fgr.p);
  fgr.perturb[1:fgr.zDim] = perturb*randn(fgr.zDim)
  copy!(fgr.X, x0)

  numericRootGenericRandomizedFnc!( fgr, perturb=perturb, testshuffle=testshuffle )
  fgr.Y
end


# function oneparams!(res::Array{Float64}, p::GenericWrapParam)
#   p.usrfnc!(res, p.params[1][:,p.particleidx])
# end
# function twoparams!(res::Array{Float64}, p::GenericWrapParam)
#   p.usrfnc!(res, p.params[1][:,p.particleidx], p.params[2][:,p.particleidx])
# end
# function threeparams!(res::Array{Float64}, p::GenericWrapParam)
#   p.usrfnc!(res, p.params[1][:,p.particleidx], p.params[2][:,p.particleidx], p.params[3][:,p.particleidx])
# end
# function fourparams!(res::Array{Float64}, p::GenericWrapParam)
#   p.usrfnc!(res, p.params[1][:,p.particleidx], p.params[2][:,p.particleidx], p.params[3][:,p.particleidx], p.params[4][:,p.particleidx])
# end


#
