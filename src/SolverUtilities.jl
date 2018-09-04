
function numericRoot(residFnc::Function, measurement, parameters, x0::Vector{Float64})
  return (nlsolve(   (res, X) -> residFnc(res, measurement, parameters, X), x0 )).zero
end


function freshSamples!(ccwl::CommonConvWrapper, N::Int=1)
  ccwl.measurement = getSample(ccwl.usrfnc!, N)
  nothing
end

function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
  # n = length(X)
  Y = deepcopy(Alt)
  for i in 1:d
    Y[p[i]] = X[i]
  end
  return Y
end

# Shuffle incoming X into random positions in fr.Y
# shuffled fr.Y will be placed back into fr.X[:,fr.gwp.particleidx] upon fr.gwp.usrfnc(x, res)
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


function (ccw::CommonConvWrapper)(res::Vector{Float64}, x::Vector{Float64})
  shuffleXAltD!(ccw, x)
  ccw.params[ccw.varidx][:, ccw.cpt[Threads.threadid()].particleidx] = ccw.cpt[Threads.threadid()].Y
  ret = ccw.usrfnc!(res, ccw.cpt[Threads.threadid()].factormetadata, ccw.cpt[Threads.threadid()].particleidx, ccw.measurement, ccw.params[ccw.cpt[Threads.threadid()].activehypo]...) # optmize the view here, re-use the same memory
  return ret
end


function (ccw::CommonConvWrapper)(x::Vector{Float64})
  ccw.params[ccw.varidx][:, ccw.cpt[Threads.threadid()].particleidx] = x #ccw.Y
  ccw.usrfnc!(ccw.cpt[Threads.threadid()].res, ccw.cpt[Threads.threadid()].factormetadata, ccw.cpt[Threads.threadid()].particleidx, ccw.measurement, ccw.params[ccw.cpt[Threads.threadid()].activehypo]...)
end



function numericRootGenericRandomizedFnc!(
            ccwl::CommonConvWrapper{T};
            perturb::Float64=1e-10,
            testshuffle::Bool=false ) where {T <: FunctorPairwiseMinimize}
  #
  fill!(ccwl.cpt[Threads.threadid()].res, 0.0) # 1:frl.xDim
  r = optimize( ccwl, ccwl.cpt[Threads.threadid()].X[:, ccwl.cpt[Threads.threadid()].particleidx] ) # ccw.gg
  # TODO -- clearly lots of optmization to be done here
  ccwl.cpt[Threads.threadid()].Y[:] = r.minimizer
  ccwl.cpt[Threads.threadid()].X[:,ccwl.cpt[Threads.threadid()].particleidx] = ccwl.cpt[Threads.threadid()].Y
  nothing
end

## TODO desperately needs cleaning up and refactoring
# Solve free variable x by root finding residual function fgr.usrfnc(x, res)
# randomly shuffle x dimensions if underconstrained by measurement z dimensions
# small random perturbation used to prevent trivial solver cases, div by 0 etc.
# result stored in fgr.Y
# fr.X must be set to memory ref the param[varidx] being solved, at creation of fr
function numericRootGenericRandomizedFnc!(
            ccwl::CommonConvWrapper{T};
            perturb::Float64=1e-10,
            testshuffle::Bool=false ) where {T <: FunctorPairwise}
  #
  # ststr = "thrid=$(Threads.threadid()), zDim=$(ccwl.zDim), xDim=$(ccwl.xDim)\n"
  # ccall(:jl_, Void, (Any,), ststr)
  if ccwl.zDim < ccwl.xDim && !ccwl.partial || testshuffle
    # less measurement dimensions than variable dimensions -- i.e. shuffle
    shuffle!(ccwl.cpt[Threads.threadid()].p)
    for i in 1:ccwl.xDim
      ccwl.cpt[Threads.threadid()].perturb[1:ccwl.zDim] = perturb*randn(ccwl.zDim)
      ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx] += ccwl.cpt[Threads.threadid()].perturb
      r = nlsolve(  ccwl,
                    ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx] # this is x0
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
    # ccall(:jl_, Void, (Any,), str)
    r = nlsolve( ccwl, ccwl.cpt[Threads.threadid()].X[1:ccwl.xDim,ccwl.cpt[Threads.threadid()].particleidx] )
    # ccall(:jl_, Void, (Any,), "nlsolve.zero=$(r.zero)")
    if sum(isnan.(( r ).zero)) == 0
      ccwl.cpt[Threads.threadid()].Y[1:ccwl.xDim] = ( r ).zero
    else
      # TODO print this output as needed
      # warn("got NaN, ccwl.cpt[Threads.threadid()].particleidx = $(ccwl.cpt[Threads.threadid()].particleidx), r=$(r)")
      str = "ccw.thrid_=$(Threads.threadid()), got NaN, ccwl.cpt[Threads.threadid()].particleidx = $(ccwl.cpt[Threads.threadid()].particleidx), r=$(r)\n"
      ccall(:jl_, Void, (Any,), str)
      ccall(:jl_, Void, (Any,), ccwl.usrfnc!)
      for thatlen in 1:length(ccwl.params)
        str = "thatlen=$thatlen, ccwl.params[thatlen][:, ccwl.cpt[Threads.threadid()].particleidx]=$(ccwl.params[thatlen][:, ccwl.cpt[Threads.threadid()].particleidx])\n"
        ccall(:jl_, Void, (Any,), str)
      end
    end
  elseif ccwl.partial
    # improve memory management in this function
    ccwl.cpt[Threads.threadid()].p[1:length(ccwl.usrfnc!.partial)] = Int[ccwl.usrfnc!.partial...] # TODO -- move this line up and out of inner loop
    r = nlsolve(  ccwl,
                  ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx] # this is x0
               )
    shuffleXAltD!( ccwl, r.zero )
  else
    ccall(:jl_, Void, (Any,), "ERROR: Unresolved numeric solve case")
    error("Unresolved numeric solve case")
  end
  ccwl.cpt[Threads.threadid()].X[:,ccwl.cpt[Threads.threadid()].particleidx] = ccwl.cpt[Threads.threadid()].Y
  nothing
end


"""
    $(SIGNATURES)

Perform multimodal incremental smoothing and mapping (mm-iSAM) computations over given factor graph `fgl::FactorGraph` on the local computer.  A pdf of the Bayes (Junction) tree will be generated in the working folder with `drawpdf=true`
"""
function batchSolve!(fgl::FactorGraph; drawpdf::Bool=false)
  if fgl.isfixedlag
      println("Quasi fixed-lag is enabled (a feature currently in testing)!")
      fifoFreeze!(fgl)
  end
  tree = wipeBuildNewTree!(fgl, drawpdf=drawpdf)
  inferOverTree!(fgl, tree)
  tree
end

"""
    $(SIGNATURES)

Update the frozen node
"""
function setfreeze!(fgl::FactorGraph, sym::Symbol)
  if !isInitialized(fgl, sym)
    warn("Vertex $(sym) is not initialized, and won't be frozen at this time.")
    return nothing
  end
  vert = getVert(fgl, sym)
  data = getData(vert)
  data.isfrozen = true

  nothing
end


"""
    $(SIGNATURES)

Freeze nodes that are older than the quasi fixed-lag length defined by `fg.qfl`, according to `fg.fifo` ordering.

Future:
- Allow different freezing strategies beyond fifo.
"""
function fifoFreeze!(fgl::FactorGraph)
  if fgl.qfl == 0
    warn("Quasi fixed-lag is enabled buyt QFL horizon is zero. Please set a valid window with FactoGraph.qfl")
  end

  tofreeze = fgl.fifo[1:(end-fgl.qfl)]
  if length(tofreeze) == 0
      info("[fifoFreeze] QFL - no nodes to freeze.")
      return nothing
  end
  info("[fifoFreeze] QFL - Freezing nodes $(tofreeze[1]) -> $(tofreeze[end]).")
  setfreeze!.(fgl, tofreeze)
  nothing
end


#
