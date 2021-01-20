

export numericSolutionCCW!

# internal function to dispatch view on either vector or matrix, rows are dims and samples are columns
_viewdim1or2(other, ind1, ind2) = other
_viewdim1or2(arr::AbstractVector, ind1, ind2) = view(arr, ind2)
_viewdim1or2(arr::AbstractMatrix, ind1, ind2) = view(arr, ind1, ind2)


function _buildCalcFactorMixture( ccwl::CommonConvWrapper,
                                  _fmd_,
                                  smpid,
                                  measurement_,
                                  varParams )
  #
  CalcFactor( ccwl.usrfnc!, _fmd_, smpid, 
              length(measurement_), measurement_, varParams)
end


function _buildCalcFactorMixture( ccwl::CommonConvWrapper{Mixture{N_,F,S,T}},
                                  _fmd_,
                                  smpid,
                                  measurement_,
                                  varParams ) where {N_,F <: FunctorInferenceType,S,T}
  #
  # just a passthrough similar to pre-v0.20
  CalcFactor( ccwl.usrfnc!.mechanics, _fmd_, smpid, 
              length(measurement_), measurement_, varParams)
end

"""
    $SIGNATURES
Internal function to build lambda pre-objective function for finding factor residuals. 

Notes  
- Unless passed in as separate arguments, this assumes already valid:
  - `cpt_.p`
  - `cpt_.activehypo`
  - `cpt_.factormetadata`
  - `ccwl.params`
  - `ccwl.measurement`

DevNotes
- TODO refactor relationship and common fields between (CCW, FMd, CPT, CalcFactor)
"""
function _buildCalcFactorLambdaSample(ccwl::CommonConvWrapper,
                                      smpid::Int,
                                      cpt_::ConvPerThread = ccwl.cpt[Threads.threadid()],
                                      target::AbstractVector = view(ccwl.params[ccwl.varidx], cpt_.p, smpid),
                                      measurement_ = ccwl.measurement,
                                      fmd_::FactorMetadata = cpt_.factormetadata  )
  #

  # build a view to the decision variable memory
  varParams = view(ccwl.params, cpt_.activehypo)
  
  # prepare fmd according to hypo selection
  # FIXME must refactor (memory waste)
  _fmd_ = FactorMetadata( view(fmd_.fullvariables, cpt_.activehypo), 
                          view(fmd_.variablelist, cpt_.activehypo),
                          varParams, # view(fmd_.arrRef, cpt_.activehypo),
                          fmd_.solvefor,
                          fmd_.cachedata  )
  #
  # get the operational CalcFactor object
  cf = _buildCalcFactorMixture(ccwl, _fmd_, smpid, measurement_, varParams)
  # new dev work on CalcFactor
  # cf = CalcFactor(ccwl.usrfnc!, _fmd_, smpid, 
  #                 length(measurement_), measurement_, varParams)
  #

  # reset the residual vector
  fill!(cpt_.res, 0.0) # 1:frl.xDim

  # build static lambda
  unrollHypo! = (res) -> cf( res, (_viewdim1or2.(measurement_, :, smpid))..., (view.(varParams, :, smpid))... )

  return unrollHypo!, target
end


# internal use only, and selected out from approxDeconv functions
_solveLambdaNumericDeconv(fcttype::AbstractPrior,
                    objResX::Function,
                    residual::AbstractVector{<:Real},
                    u0::AbstractVector{<:Real}; islen1::Bool=false ) = u0
#

function _solveLambdaNumericDeconv( fcttype::AbstractRelativeRoots,
                              objResX::Function,
                              residual::AbstractVector{<:Real},
                              u0::AbstractVector{<:Real};
                              islen1::Bool=false )
  #
  r = nlsolve(objResX, u0)

  #
  return r.zero
end

function _solveLambdaNumericDeconv( fcttype::AbstractRelativeMinimize,
                              objResX::Function,
                              residual::AbstractVector{<:Real},
                              u0::AbstractVector{<:Real};
                              islen1::Bool=false )
                              # retries::Int=3 )
  #
  r = if islen1
    optimize((x) -> objResX(residual, x), u0, BFGS() )
  else
    optimize((x) -> objResX(residual, x), u0)
  end

  # 
  return r.minimizer
end



"""
    $(SIGNATURES)

Solve free variable x by root finding residual function `fgr.usrfnc(res, x)`

ccw.X must be set to memory ref the param[varidx] being solved, at creation of ccw

Notes
- Assumes `cpt_.p` is already set to desired X decision variable dimensions and size. 
- Assumes only `ccw.particleidx` will be solved for
- small random (off-manifold) perturbation used to prevent trivial solver cases, div by 0 etc.
- Also incorporates the active hypo lookup

DevNotes
- TODO testshuffle is now obsolete, should be removed
- TODO perhaps consolidate perturbation with inflation or nullhypo
"""
function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},
                                          CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  ) where {N_,F<:AbstractRelativeMinimize,S,T}
  #
  thrid = Threads.threadid()
  smpid = ccwl.cpt[thrid].particleidx
  cpt_ = ccwl.cpt[thrid]
  
  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, smpid, cpt_)
  
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
  
  # Check for NaNs
  if sum(isnan.(( r ).minimizer)) != 0
    @error "$(ccwl.usrfnc!), ccw.thrid_=$(thrid), got NaN, smpid = $(smpid), r=$(r)\n"
    return nothing
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

  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, smpid, cpt_)

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
    return nothing
  end

  # insert result back at the correct variable element location
  cpt_.X[:,smpid] .= ( r ).zero

  nothing
end







#