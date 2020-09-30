

function _numericSolutionCCW!(fnc::InstanceType{AbstractRelativeMinimize}, 
                              ccwl::CommonConvWrapper;
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  )
  #

  # @show fnc

  fill!(ccwl.cpt[Threads.threadid()].res, 0.0) # 1:frl.xDim

  # r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx] ) # ccw.gg
  # TODO -- clearly lots of optmization to be done here
  if !ccwl.partial
    r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx] ) # ccw.gg
    ccwl.cpt[Threads.threadid()].Y .= r.minimizer
    ccwl.cpt[Threads.threadid()].X[:,ccwl.cpt[Threads.threadid()].particleidx] .= ccwl.cpt[Threads.threadid()].Y
  else
    ccwl.cpt[Threads.threadid()].p = Int[ccwl.usrfnc!.partial...]
    # ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx]
    r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[ ccwl.cpt[Threads.threadid()].p,
                                                        ccwl.cpt[Threads.threadid()].particleidx] )
    #
    ccwl.cpt[Threads.threadid()].Y = r.minimizer
    ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p,ccwl.cpt[Threads.threadid()].particleidx] .= ccwl.cpt[Threads.threadid()].Y
  end
  
  nothing
end



function _numericSolutionCCW!(fnc::InstanceType{AbstractRelativeRoots}, 
                              ccwl::CommonConvWrapper;
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  )
  #

  ## TODO desperately needs cleaning up and refactoring
  # ststr = "thrid=$(Threads.threadid()), zDim=$(ccwl.zDim), xDim=$(ccwl.xDim)\n"
  # ccall(:jl_, Nothing, (Any,), ststr)
  if ccwl.zDim < ccwl.xDim && !ccwl.partial || testshuffle
    # less measurement dimensions than variable dimensions -- i.e. shuffle
    shuffle!(ccwl.cpt[Threads.threadid()].p)
    for i in 1:ccwl.xDim
      ccwl.cpt[Threads.threadid()].perturb[1:ccwl.zDim] = perturb*randn(ccwl.zDim)
      ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx] += ccwl.cpt[Threads.threadid()].perturb
      r = nlsolve(  ccwl,
                    ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx], # this is x0
                    inplace=true
                  )
      if r.f_converged
        shuffleXAltD!( ccwl, r.zero )
        break;
      else
        # TODO -- report on this bottleneck, useful for optimization of code
        # @show i, ccwl.p, ccwl.xDim, ccwl.zDim
        temp = ccwl.cpt[Threads.threadid()].p[end]
        ccwl.cpt[Threads.threadid()].p[2:end] = ccwl.cpt[Threads.threadid()].p[1:(end-1)]
        ccwl.cpt[Threads.threadid()].p[1] = temp
        if i == ccwl.xDim
          error("numericRootGenericRandomizedFnc could not converge, i=$(i), ccwl.usrfnc!=$(typeof(ccwl.usrfnc!))")
        end
      end
    end
    #shuffleXAltD!( ccwl, r.zero ) # moved up
  elseif ccwl.zDim >= ccwl.xDim && !ccwl.partial
    # equal or more measurement dimensions than variable dimensions -- i.e. don't shuffle
    ccwl.cpt[Threads.threadid()].perturb[1:ccwl.xDim] = perturb*randn(ccwl.xDim)
    ccwl.cpt[Threads.threadid()].X[1:ccwl.xDim, ccwl.cpt[Threads.threadid()].particleidx] += ccwl.cpt[Threads.threadid()].perturb[1:ccwl.xDim] # moved up

    # str = "nlsolve, thrid_=$(Threads.threadid()), partidx=$(ccwl.cpt[Threads.threadid()].particleidx), X=$(ccwl.cpt[Threads.threadid()].X[1:ccwl.xDim,ccwl.cpt[Threads.threadid()].particleidx])"
    # ccall(:jl_, Nothing, (Any,), str)
    r = nlsolve( ccwl, ccwl.cpt[Threads.threadid()].X[1:ccwl.xDim,ccwl.cpt[Threads.threadid()].particleidx], inplace=true )
    # ccall(:jl_, Nothing, (Any,), "nlsolve.zero=$(r.zero)")
    if sum(isnan.(( r ).zero)) == 0
      ccwl.cpt[Threads.threadid()].Y[1:ccwl.xDim] = ( r ).zero
    else
      # TODO print this output as needed
      str = "ccw.thrid_=$(Threads.threadid()), got NaN, ccwl.cpt[Threads.threadid()].particleidx = $(ccwl.cpt[Threads.threadid()].particleidx), r=$(r)\n"
      @info str
      ccall(:jl_, Nothing, (Any,), str)
      ccall(:jl_, Nothing, (Any,), ccwl.usrfnc!)
      for thatlen in 1:length(ccwl.params)
        str = "thatlen=$thatlen, ccwl.params[thatlen][:, ccwl.cpt[Threads.threadid()].particleidx]=$(ccwl.params[thatlen][:, ccwl.cpt[Threads.threadid()].particleidx])\n"
        ccall(:jl_, Nothing, (Any,), str)
        @warn str
      end
    end
  elseif ccwl.partial
    # improve memory management in this function
    # TODO -- move this line up and out of inner loop
    ccwl.cpt[Threads.threadid()].p = Int[ccwl.usrfnc!.partial...] # [1:length(ccwl.usrfnc!.partial)]
    r = nlsolve(  ccwl,
                  ccwl.cpt[Threads.threadid()].X[ ccwl.cpt[Threads.threadid()].p,
                                                  ccwl.cpt[Threads.threadid()].particleidx],  # p[1:ccwl.zDim]# this is x0
                  inplace=true
                )
    shuffleXAltD!( ccwl, r.zero )
  else
    ccall(:jl_, Nothing, (Any,), "ERROR: Unresolved numeric solve case")
    error("Unresolved numeric solve case")
  end
  ccwl.cpt[Threads.threadid()].X[:,ccwl.cpt[Threads.threadid()].particleidx] = ccwl.cpt[Threads.threadid()].Y

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
function numericRootGenericRandomizedFnc!(ccwl::CommonConvWrapper{MixtureRelative{N,F,S,T}};
                                          perturb::Float64=1e-10,
                                          testshuffle::Bool=false ) where 
                                              {N,F<:AbstractRelative,S,T <: Tuple}
  #
  _numericSolutionCCW!(F, ccwl,perturb=perturb, testshuffle=testshuffle)
end


function numericRootGenericRandomizedFnc!(ccwl::CommonConvWrapper{T};
                                          perturb::Float64=1e-10,
                                          testshuffle::Bool=false ) where 
                                              {T <: AbstractRelative}
  #
  _numericSolutionCCW!(T, ccwl, perturb=perturb, testshuffle=testshuffle)
end






#