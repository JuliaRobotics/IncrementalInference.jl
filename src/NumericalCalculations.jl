

export numericSolutionCCW!

# internal function to dispatch view on either vector or matrix, rows are dims and samples are columns
_viewdim1or2(other, ind1, ind2) = other
_viewdim1or2(arr::AbstractVector, ind1, ind2) = view(arr, ind2)
_viewdim1or2(arr::AbstractMatrix, ind1, ind2) = view(arr, ind1, ind2)


function _buildCalcFactorLambdaSample(ccwl::CommonConvWrapper, 
                                      smpid::Int,
                                      thrid::Int=Threads.threadid() )
  #
  # assumes already valid `ccwl.cpt[thrid].p`
  # build a view to the decision variable memory
  varParams = view(ccwl.params, ccwl.cpt[thrid].activehypo)
  target = view(ccwl.params[ccwl.varidx], ccwl.cpt[thrid].p, smpid)
  
  # prepare fmd according to hypo selection
  # FIXME must refactor (memory waste)
  fmd = ccwl.cpt[thrid].factormetadata
  fmd_ = FactorMetadata(view(fmd.fullvariables, ccwl.cpt[thrid].activehypo), 
                        view(fmd.variablelist, ccwl.cpt[thrid].activehypo),
                        view(fmd.arrRef, ccwl.cpt[thrid].activehypo), # FIXME arrRef is likely duplicate of varParams
                        fmd.solvefor,
                        fmd.cachedata  )
  #
  # new dev work on CalcFactor
  cf = CalcFactor(ccwl.usrfnc!, fmd_, smpid, 
                  length(ccwl.measurement), ccwl.measurement, varParams)
  #

  # reset the residual vector
  fill!(ccwl.cpt[thrid].res, 0.0) # 1:frl.xDim

  # build static lambda
  unrollHypo! = (res) -> cf( res, (_viewdim1or2.(ccwl.measurement, :, smpid))..., (view.(varParams, :, smpid))... )

  return unrollHypo!, target
end



function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  )where {N_,F<:AbstractRelativeMinimize,S,T}
  #
  thrid = Threads.threadid()
  smpid = ccwl.cpt[thrid].particleidx
  
  # FIXME, can/should do this at the creation of CPT
  # indices should be permuted for Minimize
  # which elements of the variable dimension should be used as decision variables
  ccwl.cpt[thrid].p = Int[ (ccwl.partial ? ccwl.usrfnc!.partial : 1:ccwl.xDim)... ]
  
  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, smpid, thrid)
    # # assumes already valid `ccwl.cpt[thrid].p`
    # # build a view to the decision variable memory
    # varParams = view(ccwl.params, ccwl.cpt[thrid].activehypo)
    # target = view(ccwl.params[ccwl.varidx], ccwl.cpt[thrid].p, smpid)
    
    # # prepare fmd according to hypo selection
    # # FIXME must refactor (memory waste)
    # fmd = ccwl.cpt[thrid].factormetadata
    # fmd_ = FactorMetadata(view(fmd.fullvariables, ccwl.cpt[thrid].activehypo), 
    # view(fmd.variablelist, ccwl.cpt[thrid].activehypo),
    #                         view(fmd.arrRef, ccwl.cpt[thrid].activehypo), # FIXME arrRef is likely duplicate of varParams
    #                         fmd.solvefor,
    #                         fmd.cachedata  )
    # #
    # # new dev work on CalcFactor
    # cf = CalcFactor(ccwl.usrfnc!, fmd_, smpid, 
    # length(ccwl.measurement), ccwl.measurement, varParams)
    # #

    # # reset the residual vector
    # fill!(ccwl.cpt[thrid].res, 0.0) # 1:frl.xDim

    # # build static lambda
    # unrollHypo = (res) -> cf( res, (_viewdim1or2.(ccwl.measurement, :, smpid))..., (view.(varParams, :, smpid))... )
    # # unrollHypo = () -> ccwl.usrfnc!(ccwl.cpt[thrid].res,fmd_,smpid,ccwl.measurement,ccwl.params[ccwl.cpt[thrid].activehypo]...)
  
  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  _hypoObj = (x) -> (target.=x; unrollHypo!(ccwl.cpt[thrid].res) )
  # _hypoObj = (x) -> (target.=x; cf( ccwl.cpt[thrid].res ) )
  
  # cannot Nelder-Mead on 1dim
  islen1 = length(ccwl.cpt[thrid].X[:, smpid]) == 1 || ccwl.partial
  # do the parameter search over defined decision variables using Minimization
  r = if islen1
    # init value must also be permuted according to .p
    Optim.optimize( _hypoObj, ccwl.cpt[thrid].X[ccwl.cpt[thrid].p, smpid], BFGS() )
  else
    Optim.optimize( _hypoObj, ccwl.cpt[thrid].X[ccwl.cpt[thrid].p, smpid] )
  end
  
  # insert result back at the correct variable element location
  ccwl.cpt[thrid].X[ccwl.cpt[thrid].p,smpid] .= r.minimizer
  
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

  # FIXME, can/should do this at the creation of CPT
  # indices should NOT be permuted for Roots
  # which elements of the variable dimension should be used as decision variables
  ccwl.cpt[thrid].p = Int[ 1:ccwl.xDim; ] # change to `:`, type stability concern

  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, smpid, thrid)
    # # build a view to the decision variable memory
    # varParams = view(ccwl.params, ccwl.cpt[thrid].activehypo)
    # target = view(ccwl.params[ccwl.varidx], ccwl.cpt[thrid].p, smpid )
    
    # # prepare fmd according to hypo selection
    # # FIXME must refactor (memory waste)
    # fmd = ccwl.cpt[thrid].factormetadata
    # fmd_ = FactorMetadata(view(fmd.fullvariables, ccwl.cpt[thrid].activehypo), 
    #                       view(fmd.variablelist, ccwl.cpt[thrid].activehypo),
    #                       view(fmd.arrRef, ccwl.cpt[thrid].activehypo), # FIXME likely duplicate of varParams
    #                       fmd.solvefor,
    #                       fmd.cachedata  )
    # #
    # # new dev work on CalcFactor
    # cf = CalcFactor(ccwl.usrfnc!, fmd_, smpid, 
    #                 length(ccwl.measurement), ccwl.measurement, varParams)

    # # build static lambda
    # unrollHypo! = (res) -> cf( res, (_viewdim1or2.(ccwl.measurement, :, smpid))..., (view.(varParams, :, smpid))... )


  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  _hypoObj = (res,x) -> (target.=x; unrollHypo!(res))
  # unrollHypo! = (res) -> cf( res )
  # unrollHypo! = (res) -> ccwl.usrfnc!(res, fmd_, smpid, ccwl.measurement, ccwl.params[ccwl.cpt[thrid].activehypo]...)

  # NOTE small off-manifold perturbation is a numerical workaround only
  # use all element dimensions : ==> 1:ccwl.xDim
  target .+= perturb*randn(length(target))
  # ccwl.cpt[thrid].perturb[:] = perturb*randn(ccwl.xDim)
  # ccwl.cpt[thrid].X[:, smpid] += ccwl.cpt[thrid].perturb

  # do the parameter search over defined decision variables using Root finding
  r = NLsolve.nlsolve( _hypoObj, ccwl.cpt[thrid].X[:,smpid], inplace=true )
  
  # Check for NaNs
  if sum(isnan.(( r ).zero)) != 0
    @info "ccw.thrid_=$(thrid), got NaN, smpid = $(smpid), r=$(r)\n"
    for thatlen in 1:length(ccwl.params)
      @warn "thatlen=$thatlen, ccwl.params[thatlen][:, smpid]=$(ccwl.params[thatlen][:, smpid])\n"
    end
  end

  # insert result back at the correct variable element location
  ccwl.cpt[thrid].X[:,smpid] = ( r ).zero

  nothing
end







#