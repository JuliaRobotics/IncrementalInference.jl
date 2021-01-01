

export numericSolutionCCW!

# internal function to dispatch view on either vector or matrix, rows are dims and samples are columns
_viewdim1or2(other, ind1, ind2) = other
_viewdim1or2(arr::AbstractVector, ind1, ind2) = view(arr, ind2)
_viewdim1or2(arr::AbstractMatrix, ind1, ind2) = view(arr, ind1, ind2)

"""
    $SIGNATURES
Internal function to build lambda pre-objective function for finding factor residuals. 

Notes  
- assumes already valid `cpt_.p`

DevNotes
- TODO refactor relationship and common fields between (CCW, FMd, CPT, CalcFactor)
"""
function _buildCalcFactorLambdaSample(ccwl::CommonConvWrapper,
                                      cpt_::ConvPerThread,
                                      smpid::Int  )
  #

  # build a view to the decision variable memory
  varParams = view(ccwl.params, cpt_.activehypo)
  target = view(ccwl.params[ccwl.varidx], cpt_.p, smpid)
  
  # prepare fmd according to hypo selection
  # FIXME must refactor (memory waste)
  fmd = cpt_.factormetadata
  fmd_ = FactorMetadata(view(fmd.fullvariables, cpt_.activehypo), 
                        view(fmd.variablelist, cpt_.activehypo),
                        view(fmd.arrRef, cpt_.activehypo), # FIXME arrRef is likely duplicate of varParams
                        fmd.solvefor,
                        fmd.cachedata  )
  #
  # new dev work on CalcFactor
  cf = CalcFactor(ccwl.usrfnc!, fmd_, smpid, 
                  length(ccwl.measurement), ccwl.measurement, varParams)
  #

  # reset the residual vector
  fill!(cpt_.res, 0.0) # 1:frl.xDim

  # build static lambda
  unrollHypo! = (res) -> cf( res, (_viewdim1or2.(ccwl.measurement, :, smpid))..., (view.(varParams, :, smpid))... )

  return unrollHypo!, target
end



function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  ) where {N_,F<:AbstractRelativeMinimize,S,T}
  #
  thrid = Threads.threadid()
  smpid = ccwl.cpt[thrid].particleidx
  cpt_ = ccwl.cpt[thrid]
  
  # FIXME, can/should do this at the creation of CPT
  # indices should be permuted for Minimize
  # which elements of the variable dimension should be used as decision variables
  # cpt_.p = Int[ (ccwl.partial ? ccwl.usrfnc!.partial : 1:ccwl.xDim)... ]
  
  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, cpt_, smpid)
  
  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  _hypoObj = (x) -> (target.=x; unrollHypo!(cpt_.res) )
  
  # cannot Nelder-Mead on 1dim
  islen1 = length(cpt_.X[:, smpid]) == 1 || ccwl.partial
  # do the parameter search over defined decision variables using Minimization
  r = if islen1
    # init value must also be permuted according to .p
    Optim.optimize( _hypoObj, cpt_.X[cpt_.p, smpid], BFGS() )
  else
    Optim.optimize( _hypoObj, cpt_.X[cpt_.p, smpid] )
  end
  
  # insert result back at the correct variable element location
  cpt_.X[cpt_.p,smpid] .= r.minimizer
  
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

Solve free variable x by root finding residual function `fgr.usrfnc(res, x)`

ccw.X must be set to memory ref the param[varidx] being solved, at creation of ccw

Notes
- Assumes only `ccw.particleidx` will be solved for
- small random (off-manifold) perturbation used to prevent trivial solver cases, div by 0 etc.
- Also incorporates the active hypo lookup

DevNotes
- TODO testshuffle is now obsolete, should be removed
- TODO perhaps consolidate perturbation with inflation or nullhypo
"""
function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  ) where {N_,F<:AbstractRelativeRoots,S,T}
  #
  
  # FIXME, move this check higher and out of smpid loop
  if ccwl.zDim < ccwl.xDim && !ccwl.partial || testshuffle || ccwl.partial
    error("<:AbstractRelativeRoots factors with less measurement dimensions than variable dimensions have been discontinued, easy conversion to <:AbstractRelativeMinimize is the better option.")
  elseif !( ccwl.zDim >= ccwl.xDim && !ccwl.partial )
    error("Unresolved numeric <:AbstractRelativeRoots solve case")
  end

  thrid = Threads.threadid()
  smpid = ccwl.cpt[thrid].particleidx
  cpt_ = ccwl.cpt[thrid]

  # FIXME, can/should do this at the creation of CPT
  # indices should NOT be permuted for Roots
  # which elements of the variable dimension should be used as decision variables
  # TODO perhaps rename `.p` to `.decisionVarDims`
  # cpt_.p = Int[ 1:ccwl.xDim; ] # change to `:`, type stability concern

  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, cpt_, smpid)

  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  _hypoObj = (res,x) -> (target.=x; unrollHypo!(res))

  # NOTE small off-manifold perturbation is a numerical workaround only
  # use all element dimensions : ==> 1:ccwl.xDim
  target .+= perturb*randn(length(target))

  # do the parameter search over defined decision variables using Root finding
  r = NLsolve.nlsolve( _hypoObj, cpt_.X[:,smpid], inplace=true )
  
  # Check for NaNs
  if sum(isnan.(( r ).zero)) != 0
    @error "$(ccwl.usrfnc!), ccw.thrid_=$(thrid), got NaN, smpid = $(smpid), r=$(r)\n"
    # for thatlen in 1:length(ccwl.params)
    #   @warn "thatlen=$thatlen, ccwl.params[thatlen][:, smpid]=$(ccwl.params[thatlen][:, smpid])\n"
    # end
  end

  # insert result back at the correct variable element location
  cpt_.X[:,smpid] = ( r ).zero

  nothing
end







#