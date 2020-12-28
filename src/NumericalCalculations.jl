

"""
    $(SIGNATURES)

Shuffle incoming X into random positions in fr.Y.
Shuffled fr.Y will be placed back into fr.X[:,fr.gwp.particleidx] upon fr.gwp.usrfnc(x, res).

DevNotes
- TODO remove, see #1072
"""
function shuffleXAltD!(ccwl::CommonConvWrapper, X::Vector{Float64})
  # populate defaults from existing values
  for i in 1:ccwl.xDim
    ccwl.cpt[Threads.threadid()].Y[i] = ccwl.cpt[Threads.threadid()].X[i, ccwl.cpt[Threads.threadid()].particleidx]
  end
  # populate as many measurment dimensions randomly for calculation
  for i in 1:ccwl.zDim
    ccwl.cpt[Threads.threadid()].Y[ccwl.cpt[Threads.threadid()].p[i]] = X[i]
  end
  nothing
end


function (ccw::CommonConvWrapper)(res::AbstractVector{<:Real}, x::AbstractVector{<:Number})
  # TODO unclear how long this feature will persist
  shuffleXAltD!(ccw, x)
  ccw.params[ccw.varidx][:, ccw.cpt[Threads.threadid()].particleidx] = ccw.cpt[Threads.threadid()].Y
  # evaulate the user provided residual function with constructed set of parameters
  ret = ccw.usrfnc!(res,
                    ccw.cpt[Threads.threadid()].factormetadata,
                    ccw.cpt[Threads.threadid()].particleidx,
                    ccw.measurement,
                    ccw.params[ccw.cpt[Threads.threadid()].activehypo]...) # optmize the view here, re-use the same memory
  return ret
end


function (ccw::CommonConvWrapper)(x::AbstractVector{<:Number} )
  #
  # AbstractRelativeMinimize case
  # set internal memory to that of external caller value `x`, special care if partial
  if !ccw.partial
    ccw.params[ccw.varidx][:, ccw.cpt[Threads.threadid()].particleidx] .= x #ccw.Y
  else
    ccw.params[ccw.varidx][ccw.cpt[Threads.threadid()].p, ccw.cpt[Threads.threadid()].particleidx] .= x #ccw.Y
  end
  # evaluate the user provided residual function with constructed set of parameters
  ccw.usrfnc!(ccw.cpt[Threads.threadid()].res,
              ccw.cpt[Threads.threadid()].factormetadata,
              ccw.cpt[Threads.threadid()].particleidx,
              ccw.measurement,
              ccw.params[ccw.cpt[Threads.threadid()].activehypo]...)
end

# should only be calling a new arg list according to activehypo at start of particle
# Try calling an existing lambda
# sensitive to which hypo of course , see #1024
# need to shuffle content inside .cpt.fmd as well as .params accordingly
#
#
#
#
#

# function (cf::CalcResidual{<:LinearRelative})(res::Vector, z, xi, xj)
#   # cf.metadata.variablelist...
#   # cf.metadata.targetvariable
#   # cf.metadata.usercache
#   # generic on-manifold residual function 
#   res .= distance(z, distance(xj, xi))
# end

function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  )where {N_,F<:AbstractRelativeMinimize,S,T}
  #
  # fnc::InstanceType{AbstractRelativeMinimize}, 

  # @show fnc

  fill!(ccwl.cpt[Threads.threadid()].res, 0.0) # 1:frl.xDim

  # r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx] ) # ccw.gg
  # TODO -- clearly lots of optmization to be done here
  islen1 = length(ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx]) == 1
  moreargs = islen1 ? (BFGS(),) : ()
  if !ccwl.partial
    r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx], moreargs... )
    ccwl.cpt[Threads.threadid()].Y .= r.minimizer
    ccwl.cpt[Threads.threadid()].X[:,ccwl.cpt[Threads.threadid()].particleidx] .= ccwl.cpt[Threads.threadid()].Y
  else
    ccwl.cpt[Threads.threadid()].p = Int[ccwl.usrfnc!.partial...]
    # ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx]
    r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[ ccwl.cpt[Threads.threadid()].p, ccwl.cpt[Threads.threadid()].particleidx], BFGS() )
    #
    ccwl.cpt[Threads.threadid()].Y = r.minimizer
    ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p,ccwl.cpt[Threads.threadid()].particleidx] .= ccwl.cpt[Threads.threadid()].Y
  end
  
  nothing
end


"""
    $(SIGNATURES)

Solve free variable x by root finding residual function fgr.usrfnc(x, res)
randomly shuffle x dimensions if underconstrained by measurement z dimensions
small random perturbation used to prevent trivial solver cases, div by 0 etc.
result stored in fgr.Y
ccw.X must be set to memory ref the param[varidx] being solved, at creation of ccw

Notes
- Assumes only `ccw.particleidx` will be solved for
"""
function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  ) where {N_,F<:AbstractRelativeRoots,S,T}
  #
  thrid = Threads.threadid()

  ## TODO needs cleaning up and refactoring
  if ccwl.zDim >= ccwl.xDim && !ccwl.partial
    # equal or more measurement dimensions than variable dimensions -- i.e. don't shuffle
    ccwl.cpt[thrid].perturb[1:ccwl.xDim] = perturb*randn(ccwl.xDim)
    ccwl.cpt[thrid].X[1:ccwl.xDim, ccwl.cpt[thrid].particleidx] += ccwl.cpt[thrid].perturb[1:ccwl.xDim] # moved up
    
    # str = "nlsolve, thrid_=$(thrid), partidx=$(ccwl.cpt[thrid].particleidx), X=$(ccwl.cpt[thrid].X[1:ccwl.xDim,ccwl.cpt[thrid].particleidx])"
    r = nlsolve( ccwl, ccwl.cpt[thrid].X[1:ccwl.xDim,ccwl.cpt[thrid].particleidx], inplace=true )
    if sum(isnan.(( r ).zero)) == 0
      ccwl.cpt[thrid].Y[1:ccwl.xDim] = ( r ).zero
    else
      @info "ccw.thrid_=$(thrid), got NaN, ccwl.cpt[thrid].particleidx = $(ccwl.cpt[thrid].particleidx), r=$(r)\n"
      for thatlen in 1:length(ccwl.params)
        @warn "thatlen=$thatlen, ccwl.params[thatlen][:, ccwl.cpt[thrid].particleidx]=$(ccwl.params[thatlen][:, ccwl.cpt[thrid].particleidx])\n"
      end
    end
  elseif ccwl.zDim < ccwl.xDim && !ccwl.partial || testshuffle || ccwl.partial
    error("<:AbstractRelativeRoots factors with less measurement dimensions than variable dimensions have been discontinued, easy conversion to <:AbstractRelativeMinimize is the better option.")
  else
    error("Unresolved numeric solve case")
  end
  ccwl.cpt[thrid].X[:,ccwl.cpt[thrid].particleidx] = ccwl.cpt[thrid].Y

  nothing
end







#