
function numericRoot(residFnc::Function, measurement, parameters, x0::Vector{Float64})
  return (nlsolve(   (res, X) -> residFnc(res, measurement, parameters, X), x0 )).zero
end


function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
  # n = length(X)
  Y = deepcopy(Alt)
  for i in 1:d
    Y[p[i]] = X[i]
  end
  return Y
end


function (p::GenericWrapParam)(res, x)
  # TODO -- move to inner lambda that is defined once against p.params...
  # approximates by not considering cross indices among parameters
  # @show length(p.params), p.varidx, p.particleidx, size(x), size(res), size(p.measurement)
  p.params[p.varidx][:, p.particleidx] = x
  # p.usrfnc!(res, p.particleidx, p.measurement, p.params...)
  # who are active hypotheses?  p.params[p.activehypo]...
  p.usrfnc!(res, p.factormetadata, p.particleidx, p.measurement, view(p.params,p.activehypo)...)
  # p.usrfnc!(res, p.factormetadata, p.particleidx, p.measurement, p.params[p.activehypo]...)
end

# Shuffle incoming X into random permutation in fr.Y
# shuffled fr.Y will be placed back into fr.X[:,fr.gwp.particleidx] upon fr.gwp.usrfnc(x, res)
function shuffleXAltD!(fr::FastRootGenericWrapParam, X::Vector{Float64})
  fr.Y[1:fr.xDim] = view(fr.X, 1:fr.xDim, fr.gwp.particleidx)
  # fr.Y[1:fr.xDim] = fr.X[1:fr.xDim,fr.gwp.particleidx]
  # copy!(fr.Y, fr.X[:,fr.gwp.particleidx])
  for i in 1:fr.zDim
    fr.Y[fr.p[i]] = X[i]
  end
  nothing
end
function (fr::FastRootGenericWrapParam)( res::Vector{Float64}, x::Vector{Float64} )
  shuffleXAltD!(fr, x)
  fr.gwp( res, fr.Y )
end


## TODO desperately needs cleaning up and refactoring
# Solve free variable x by root finding residual function fgr.usrfnc(x, res)
# randomly shuffle x dimensions if underconstrained by measurement z dimensions
# small random perturbation used to prevent trivial solver cases, div by 0 etc.
# result stored in fgr.Y
# fr.X must be set to memory ref the param[varidx] being solved, at creation of fr
function numericRootGenericRandomizedFnc!(
      fr::FastRootGenericWrapParam{T};
      perturb::Float64=1e-10,
      testshuffle::Bool=false ) where {T <: FunctorInferenceType}
  #
  # info("numericRootGenericRandomizedFnc! FastRootGenericWrapParam{T}")
  # @show fr.zDim, fr.xDim, fr.gwp.partial, fr.gwp.particleidx
  if fr.zDim < fr.xDim && !fr.gwp.partial || testshuffle
    shuffle!(fr.p)
    for i in 1:fr.xDim
      fr.perturb[1:fr.zDim] = perturb*randn(fr.zDim)
      fr.X[fr.p[1:fr.zDim], fr.gwp.particleidx] += fr.perturb
      r = nlsolve(  fr,
                    fr.X[fr.p[1:fr.zDim], fr.gwp.particleidx] # this is x0
                 )
      if r.f_converged
        shuffleXAltD!( fr, r.zero )
        break;
      else
        # TODO -- report on this bottleneck, useful for optimization of code
        # @show i, fr.p, fr.xDim, fr.zDim
        temp = fr.p[end]
        fr.p[2:end] = fr.p[1:(end-1)]
        fr.p[1] = temp
        if i == fr.xDim
          error("numericRootGenericRandomizedFnc could not converge, i=$(i), fr.gwp.usrfnc!=$(typeof(fr.gwp.usrfnc!))")
        end
      end
    end
    #shuffleXAltD!( fr, r.zero ) # moved up
  elseif fr.zDim >= fr.xDim && !fr.gwp.partial
    fr.perturb[1:fr.xDim] = perturb*randn(fr.xDim)
    fr.X[1:fr.xDim, fr.gwp.particleidx] += fr.perturb[1:fr.xDim] # moved up
    r = nlsolve( fr.gwp, fr.X[1:fr.xDim,fr.gwp.particleidx] )
    if sum(isnan.(( r ).zero)) == 0
      fr.Y[1:fr.xDim] = ( r ).zero
    else
      warn("got NaN, fr.gwp.particleidx = $(fr.gwp.particleidx), r=$(r)")
      @show fr.gwp.usrfnc!
      for thatlen in 1:length(fr.gwp.params)
        @show thatlen, fr.gwp.params[thatlen][:, fr.gwp.particleidx]
      end
    end
  elseif fr.gwp.partial
    # improve memory management in this function
    fr.p[1:length(fr.gwp.usrfnc!.partial)] = Int[fr.gwp.usrfnc!.partial...] # TODO -- move this line up and out of inner loop
    r = nlsolve(  fr,
                  fr.X[fr.p[1:fr.zDim], fr.gwp.particleidx] # this is x0
               )
    shuffleXAltD!( fr, r.zero )
  else
    error("Unresolved numeric solve case")
  end
  fr.X[:,fr.gwp.particleidx] = fr.Y
  nothing
end


#-------------------------------------------------------------------------------

mutable struct TestSolver <: FunctorPairwise
  testfnc::Function
end


# function numericRootGenericRandomizedFnc(
#       residFnc!::Function,
#       zDim::Int,
#       xDim::Int,
#       x0::Vector{Float64};
#       perturb::Float64=1e-5,
#       testshuffle::Bool=false   )
#   #
#   # TODO -- this only start of refactoring for inplace, more to come
#   # xDim = length(x0)
#   # fgr = FastGenericRoot{typeof(residFnc!)}(xDim, zDim, residFnc!)
#   # shuffle!(fgr.p);
#   # fgr.perturb[1:fgr.zDim] = perturb*randn(fgr.zDim)
#   # copy!(fgr.X, x0)
#
#   info("This specific `numericRootGenericRandomizedFnc(,,,)` function is a convenience function only -- do tnot use in production.")
#
#   rr = TestSolver(residFnc!)
#   gwp = GenericWrapParam{TestSolver}(rr, t, 1, 1)
#   gwp.measurement = (zeros(zDim,0), )
#   fr = FastRootGenericWrapParam{TestSolver}(gwp.params[gwp.varidx], zDim, gwp)
#
#   fr.xDim = xDim
#   gwp.particleidx = 1
#
#   numericRootGenericRandomizedFnc!( fgr, perturb=perturb, testshuffle=testshuffle )
#   fgr.Y
# end



#
