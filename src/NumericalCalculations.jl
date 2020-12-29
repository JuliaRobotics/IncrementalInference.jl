

export numericSolutionCCW!

# Also see #467 on API consolidation
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
  
  # which elements of the variable dimension should be used as decision variables
  ccwl.cpt[thrid].p = Int[ (ccwl.partial ? ccwl.usrfnc!.partial : 1:ccwl.xDim)... ]  
  # build a view to the decision variable memory
  target = view(ccwl.params[ccwl.varidx], ccwl.cpt[thrid].p, ccwl.cpt[thrid].particleidx)

  # prepare fmd according to hypo selection
  # fmd = ccwl.cpt[thrid].factormetadata
  # fmd.fullvariables = 
  # fmd.arrRef = 

  # build static lambda
  unrollHypo = () -> ccwl.usrfnc!(ccwl.cpt[thrid].res,ccwl.cpt[thrid].factormetadata,ccwl.cpt[thrid].particleidx,ccwl.measurement,ccwl.params[ccwl.cpt[thrid].activehypo]...)
  # broadcast updates original view memory location
  _hypoObj = (x) -> (target.=x; unrollHypo())
  
  # cannot Nelder-Mead on 1dim
  islen1 = length(ccwl.cpt[thrid].X[:, ccwl.cpt[thrid].particleidx]) == 1 || ccwl.partial
  # do the parameter search over defined decision variables using Minimization
  r = if islen1
    Optim.optimize( _hypoObj, ccwl.cpt[thrid].X[ ccwl.cpt[thrid].p, ccwl.cpt[thrid].particleidx], BFGS() )
  else
    Optim.optimize( _hypoObj, ccwl.cpt[thrid].X[ ccwl.cpt[thrid].p, ccwl.cpt[thrid].particleidx] )
  end
  # r = Optim.optimize( ccwl, ccwl.cpt[thrid].X[ ccwl.cpt[thrid].p, ccwl.cpt[thrid].particleidx], moreargs... )
  
  # extract the result from inference
  # ccwl.cpt[thrid].Y[:] = r.minimizer
  
  # insert result back at the correct variable element location
  ccwl.cpt[thrid].X[ccwl.cpt[thrid].p,ccwl.cpt[thrid].particleidx] .= r.minimizer #  ccwl.cpt[thrid].Y
  
  nothing
end
# brainstorming
# should only be calling a new arg list according to activehypo at start of particle
# Try calling an existing lambda
# sensitive to which hypo of course , see #1024
# need to shuffle content inside .cpt.fmd as well as .params accordingly
#



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

  # build a view to the decision variable memory
  target = view(ccwl.params[ccwl.varidx], :, ccwl.cpt[Threads.threadid()].particleidx )
  
  # build static lambda
  unrollHypo! = (res) -> ccwl.usrfnc!(res, ccwl.cpt[thrid].factormetadata, ccwl.cpt[thrid].particleidx, ccwl.measurement, ccwl.params[ccwl.cpt[thrid].activehypo]...)
  # broadcast updates original view memory location
  _hypoObj = (res,x) -> (target.=x; unrollHypo!(res))

  # do the parameter search over defined decision variables using Root finding
  r = NLsolve.nlsolve( _hypoObj, ccwl.cpt[thrid].X[:,ccwl.cpt[thrid].particleidx], inplace=true )
  # r = NLsolve.nlsolve( ccwl, ccwl.cpt[thrid].X[:,ccwl.cpt[thrid].particleidx], inplace=true )
  
  # Check for NaNs
  if sum(isnan.(( r ).zero)) == 0
    # ccwl.cpt[thrid].Y[:] = ( r ).zero
  else
    @info "ccw.thrid_=$(thrid), got NaN, ccwl.cpt[thrid].particleidx = $(ccwl.cpt[thrid].particleidx), r=$(r)\n"
    for thatlen in 1:length(ccwl.params)
      @warn "thatlen=$thatlen, ccwl.params[thatlen][:, ccwl.cpt[thrid].particleidx]=$(ccwl.params[thatlen][:, ccwl.cpt[thrid].particleidx])\n"
    end
  end

  # insert result back at the correct variable element location
  ccwl.cpt[thrid].X[:,ccwl.cpt[thrid].particleidx] = ( r ).zero # ccwl.cpt[thrid].Y

  nothing
end







#