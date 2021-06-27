

# TODO deprecate testshuffle
_checkErrorCCWNumerics(ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}}, testshuffle::Bool=false)  where {N_,F<:AbstractRelativeMinimize,S,T} = nothing

function _checkErrorCCWNumerics(ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}},
                                testshuffle::Bool=false)  where {N_,F<:AbstractRelativeRoots,S,T}
  #
  @info "ccwl zDim and xDim" ccwl.zDim ccwl.xDim
  if ccwl.zDim < ccwl.xDim && !ccwl.partial || testshuffle || ccwl.partial
    error("<:AbstractRelativeRoots factors with less measurement dimensions than variable dimensions have been discontinued, easy conversion to <:AbstractRelativeMinimize is the better option.")
  elseif !( ccwl.zDim >= ccwl.xDim && !ccwl.partial )
    error("Unresolved numeric <:AbstractRelativeRoots solve case")
  end
  nothing
end


_perturbIfNecessary(fcttype::Union{F,<:Mixture{N_,F,S,T}},
                    len::Int=1,
                    perturbation::Real=1e-10 ) where {N_,F<:AbstractRelativeMinimize,S,T} = 0
#

_perturbIfNecessary(fcttype::Union{F,<:Mixture{N_,F,S,T}},
                    len::Int=1,
                    perturbation::Real=1e-10 ) where {N_,F<:AbstractRelativeRoots,S,T} = perturbation*randn(len)
#


# internal use only, and selected out from approxDeconv functions
_solveLambdaNumeric(fcttype::AbstractPrior,
                    objResX::Function,
                    residual::AbstractVector{<:Real},
                    u0::AbstractVector{<:Real},
                    islen1::Bool=false;
                    perturb::Real=1e-10 ) = u0
#

function _solveLambdaNumeric( fcttype::Union{F,<:Mixture{N_,F,S,T}},
                              objResX::Function,
                              residual::AbstractVector{<:Real},
                              u0::AbstractVector{<:Real},
                              islen1::Bool=false )  where {N_,F<:AbstractRelativeRoots,S,T}
  #

  #
  r = NLsolve.nlsolve( (res, x) -> res .= objResX(x), u0, inplace=true) #, ftol=1e-14)

  #
  return r.zero
end

function _solveLambdaNumeric( fcttype::Union{F,<:Mixture{N_,F,S,T}},
                              objResX::Function,
                              residual::AbstractVector{<:Real},
                              u0::AbstractVector{<:Real},
                              islen1::Bool=false  )  where {N_,F<:AbstractRelativeMinimize,S,T}
                              # retries::Int=3 )
  #
  # wrt #467 allow residual to be standardize for Roots and Minimize and Parametric cases.
  r = if islen1
    Optim.optimize((x) -> (residual .= objResX(x); sum(residual.^2)), u0, Optim.BFGS() )
  else
    Optim.optimize((x) -> (residual .= objResX(x); sum(residual.^2)), u0)
  end

  # 
  return r.minimizer
end



## ================================================================================================
## Heavy dispatch for all AbstractFactor / Mixture cases below
## ================================================================================================



# internal function to dispatch view on either vector or matrix, rows are dims and samples are columns
# _viewdim1or2(other, ind...) = other
_getindextuple(tup::Tuple, ind1) = [getindex(t, ind1) for t in tup]
# _viewdim1or2(arr::AbstractMatrix, ind1, ind2) = view(arr, ind1, ind2)


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
- Unless passed in as separate arguments, this assumes already valid in `cpt_`:
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
                                      target::AbstractVector = view(ccwl.params[ccwl.varidx][smpid], cpt_.p),
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
  fill!(cpt_.res, 0.0) # Roots->xDim | Minimize->zDim

  # build static lambda
  # unrollHypo! =  ()->cf( (_viewdim1or2.(measurement_, smpid))..., (view.(varParams, smpid))... ) # :
  unrollHypo! =  ()->cf( (_getindextuple(measurement_, smpid))..., (getindex.(varParams, smpid))... ) # :

  return unrollHypo!, target
end


"""
    $(SIGNATURES)

Solve free variable x by root finding residual function `fgr.usrfnc(res, x)`.  This is the 
penultimate step before calling numerical operations to move actual estimates, which is 
done by an internally created lambda function.

ccw.X must be set to memory ref the param[varidx] being solved, at creation of ccw

Notes
- Assumes `cpt_.p` is already set to desired X decision variable dimensions and size. 
- Assumes only `ccw.particleidx` will be solved for
- small random (off-manifold) perturbation used to prevent trivial solver cases, div by 0 etc.
  - perturb is necessary for NLsolve cases, and smaller than 1e-10 will result in test failure
- Also incorporates the active hypo lookup

DevNotes
- TODO testshuffle is now obsolete, should be removed
- TODO perhaps consolidate perturbation with inflation or nullhypo
"""
function _solveCCWNumeric!( ccwl::Union{CommonConvWrapper{F},
                                        CommonConvWrapper{Mixture{N_,F,S,T}}};
                            perturb::Real=1e-10,
                            testshuffle::Bool=false  ) where {N_,F<:AbstractRelative,S,T}
  #
    # FIXME, move this check higher and out of smpid loop
  _checkErrorCCWNumerics(ccwl, testshuffle)

  #
  thrid = Threads.threadid()
  cpt_ = ccwl.cpt[thrid]
  smpid = cpt_.particleidx
  # cannot Nelder-Mead on 1dim, partial can be 1dim or more but being conservative.
  islen1 = length(cpt_.p) == 1 || ccwl.partial
  # islen1 = length(cpt_.X[:, smpid]) == 1 || ccwl.partial

  # build the pre-objective function for this sample's hypothesis selection
  unrollHypo!, target = _buildCalcFactorLambdaSample(ccwl, smpid, cpt_)
  
  # broadcast updates original view memory location
  ## using CalcFactor legacy path inside (::CalcFactor)
  _hypoObj = (x) -> (target.=x; unrollHypo!())
  
  # TODO small off-manifold perturbation is a numerical workaround only, make on-manifold requires RoME.jl #244
  # use all element dimensions : ==> 1:ccwl.xDim
  target .+= _perturbIfNecessary(getFactorType(ccwl), length(target), perturb)

  # do the parameter search over defined decision variables using Minimization
  retval = _solveLambdaNumeric(getFactorType(ccwl), _hypoObj, cpt_.res, cpt_.X[smpid][cpt_.p], islen1 )
  
  # Check for NaNs
  if sum(isnan.(retval)) != 0
    @error "$(ccwl.usrfnc!), ccw.thrid_=$(thrid), got NaN, smpid = $(smpid), r=$(retval)\n"
    return nothing
  end

  # insert result back at the correct variable element location
  cpt_.X[smpid][cpt_.p] .= retval
  
  nothing
end
# brainstorming
# should only be calling a new arg list according to activehypo at start of particle
# Try calling an existing lambda
# sensitive to which hypo of course , see #1024
# need to shuffle content inside .cpt.fmd as well as .params accordingly
#





#