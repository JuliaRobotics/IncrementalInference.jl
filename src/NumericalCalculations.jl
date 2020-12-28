

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
  # <: AbstractRelativeRoots case
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
  idxs = ccw.cpt[Threads.threadid()].p
  # set internal memory to that of external caller value `x`, special care if partial
  ccw.params[ccw.varidx][idxs, ccw.cpt[Threads.threadid()].particleidx] .= x #ccw.Y
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
  thrid = Threads.threadid()

  # reset the residual vector
  fill!(ccwl.cpt[thrid].res, 0.0) # 1:frl.xDim

  # cannot Nelder-Mead on 1dim
  islen1 = length(ccwl.cpt[thrid].X[:, ccwl.cpt[thrid].particleidx]) == 1
  moreargs = islen1 || ccwl.partial ? (BFGS(),) : () 

  # which elements of the variable dimension should be used as decision variables
  ccwl.cpt[thrid].p = Int[ (ccwl.partial ? ccwl.usrfnc!.partial : 1:ccwl.xDim)... ]

  # do the parameter search over defined decision variables using Minimization
  r = Optim.optimize( ccwl, ccwl.cpt[thrid].X[ ccwl.cpt[thrid].p, ccwl.cpt[thrid].particleidx], moreargs... )

  # extract the result from inference
  ccwl.cpt[thrid].Y = r.minimizer
  
  # insert result back in variable element
  ccwl.cpt[thrid].X[ccwl.cpt[thrid].p,ccwl.cpt[thrid].particleidx] .= ccwl.cpt[thrid].Y
  
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

DevNotes
- TODO perhaps consolidate perturbation with inflation or nullhypo
"""
function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  ) where {N_,F<:AbstractRelativeRoots,S,T}
  #
  thrid = Threads.threadid()

  if ccwl.zDim < ccwl.xDim && !ccwl.partial || testshuffle || ccwl.partial
    error("<:AbstractRelativeRoots factors with less measurement dimensions than variable dimensions have been discontinued, easy conversion to <:AbstractRelativeMinimize is the better option.")
  elseif !( ccwl.zDim >= ccwl.xDim && !ccwl.partial )
    error("Unresolved numeric <:AbstractRelativeRoots solve case")
  end
  # equal or more measurement dimensions than variable dimensions -- i.e. don't shuffle
  # if ccwl.zDim >= ccwl.xDim && !ccwl.partial
  
  # NOTE ignoring small perturbation from manifold as numerical workaround only
  # use all element dimensions : ==> 1:ccwl.xDim
  ccwl.cpt[thrid].perturb[:] = perturb*randn(ccwl.xDim)
  ccwl.cpt[thrid].X[:, ccwl.cpt[thrid].particleidx] += ccwl.cpt[thrid].perturb
  
  # do the parameter search over defined decision variables using Root finding
  r = NLsolve.nlsolve( ccwl, ccwl.cpt[thrid].X[:,ccwl.cpt[thrid].particleidx], inplace=true )
  
  # Check for NaNs
  if sum(isnan.(( r ).zero)) == 0
    ccwl.cpt[thrid].Y[:] = ( r ).zero
  else
    @info "ccw.thrid_=$(thrid), got NaN, ccwl.cpt[thrid].particleidx = $(ccwl.cpt[thrid].particleidx), r=$(r)\n"
    for thatlen in 1:length(ccwl.params)
      @warn "thatlen=$thatlen, ccwl.params[thatlen][:, ccwl.cpt[thrid].particleidx]=$(ccwl.params[thatlen][:, ccwl.cpt[thrid].particleidx])\n"
    end
  end

  # insert result back in variable element
  ccwl.cpt[thrid].X[:,ccwl.cpt[thrid].particleidx] = ccwl.cpt[thrid].Y

  nothing
end







#