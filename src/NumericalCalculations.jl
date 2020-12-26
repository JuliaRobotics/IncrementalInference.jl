

function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
  # n = length(X)
  Y = deepcopy(Alt)
  for i in 1:d
    Y[p[i]] = X[i]
  end
  return Y
end

"""
    $(SIGNATURES)

Shuffle incoming X into random positions in fr.Y.
Shuffled fr.Y will be placed back into fr.X[:,fr.gwp.particleidx] upon fr.gwp.usrfnc(x, res).
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


function (ccw::CommonConvWrapper)(res::AbstractVector{<:Real}, x::AbstractVector{<:Real})
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


function (ccw::CommonConvWrapper)(x::Vector{Float64})
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



function numericSolutionCCW!( ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}};
                              perturb::Float64=1e-10,
                              testshuffle::Bool=false  )where {N_,F<:AbstractRelativeMinimize,S,T}
  #
  # fnc::InstanceType{AbstractRelativeMinimize}, 

  # @show fnc

  fill!(ccwl.cpt[Threads.threadid()].res, 0.0) # 1:frl.xDim

  # r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx] ) # ccw.gg
  # TODO -- clearly lots of optmization to be done here
  if !ccwl.partial
    r = if length(ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx]) == 1
      optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx], BFGS() )
    else
      optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx] )
    end
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
          error("numericSolutionCCW! could not converge, i=$(i), ccwl.usrfnc!=$(typeof(ccwl.usrfnc!))")
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







#