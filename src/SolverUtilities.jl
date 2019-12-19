function fastnorm(u)
  # dest[1] = ...
  n = length(u)
  T = eltype(u)
  s = zero(T)
  @fastmath @inbounds @simd for i in 1:n
      s += u[i]^2
  end
  @fastmath @inbounds return sqrt(s)
end

function numericRoot(residFnc::Function, measurement, parameters, x0::Vector{Float64})
  # function is being deprecated
  return (nlsolve(   (res, X) -> residFnc(res, measurement, parameters, X), x0, inplace=true )).zero
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

"""
    $(SIGNATURES)

Solve free variable x by root finding residual function fgr.usrfnc(x, res)
randomly shuffle x dimensions if underconstrained by measurement z dimensions
small random perturbation used to prevent trivial solver cases, div by 0 etc.
result stored in fgr.Y
ccw.X must be set to memory ref the param[varidx] being solved, at creation of ccw
"""
function numericRootGenericRandomizedFnc!(
            ccwl::CommonConvWrapper{T};
            perturb::Float64=1e-10,
            testshuffle::Bool=false ) where {T <: FunctorPairwise}
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
    ccwl.cpt[Threads.threadid()].p[1:length(ccwl.usrfnc!.partial)] = Int[ccwl.usrfnc!.partial...] # TODO -- move this line up and out of inner loop
    r = nlsolve(  ccwl,
                  ccwl.cpt[Threads.threadid()].X[ccwl.cpt[Threads.threadid()].p[1:ccwl.zDim], ccwl.cpt[Threads.threadid()].particleidx], # this is x0
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

Set variable(s) `sym` of factor graph to be marginalized -- i.e. not be updated by inference computation.
"""
function setfreeze!(dfg::AbstractDFG, sym::Symbol)
  if !isInitialized(dfg, sym)
    @warn "Vertex $(sym) is not initialized, and won't be frozen at this time."
    return nothing
  end
  vert = DFG.getVariable(dfg, sym)
  data = solverData(vert)
  data.ismargin = true
  nothing
end
function setfreeze!(dfg::AbstractDFG, syms::Vector{Symbol})
  for sym in syms
    setfreeze!(dfg, sym)
  end
end

"""
    $(SIGNATURES)

Freeze nodes that are older than the quasi fixed-lag length defined by `fg.qfl`, according to `fg.fifo` ordering.

Future:
- Allow different freezing strategies beyond fifo.
"""
function fifoFreeze!(dfg::G)::Nothing where G <: AbstractDFG
  if DFG.getSolverParams(dfg).qfl == 0
    @warn "Quasi fixed-lag is enabled but QFL horizon is zero. Please set a valid window with FactoGraph.qfl"
  end

  # the fifo history
  tofreeze = DFG.getAddHistory(dfg)[1:(end-DFG.getSolverParams(dfg).qfl)]
  if length(tofreeze) == 0
      @info "[fifoFreeze] QFL - no nodes to freeze."
      return nothing
  end
  @info "[fifoFreeze] QFL - Freezing nodes $(tofreeze[1]) -> $(tofreeze[end])."
  setfreeze!(dfg, tofreeze)
  nothing
end

"""
    $(SIGNATURES)

Return all factors currently registered in the workspace.
"""
function getCurrentWorkspaceFactors()::Vector{Type}
    return [
        subtypes(IncrementalInference.FunctorSingleton)...,
        subtypes(IncrementalInference.FunctorPairwise)...,
        subtypes(IncrementalInference.FunctorPairwiseMinimize)...];
end

"""
    $(SIGNATURES)

Return all variables currently registered in the workspace.
"""
function getCurrentWorkspaceVariables()::Vector{Type}
    return subtypes(IncrementalInference.InferenceVariable);
end


"""
    $SIGNATURES

Reset the state of all variables in a clique to not initialized.

Notes
- resets numberical values to zeros.

Dev Notes
- TODO not all kde manifolds will initialize to zero.
"""
function resetCliqSolve!(dfg::G,
                         treel::BayesTree,
                         cliq::Graphs.ExVertex;
                         solveKey::Symbol=:default)::Nothing where G <: AbstractDFG
  #
  cda = getData(cliq)
  vars = getCliqVarIdsAll(cliq)
  for varis in vars
    resetVariable!(dfg, varis, solveKey=solveKey)
  end
  prnt = getParent(treel, cliq)
  if length(prnt) > 0
    setCliqUpInitMsgs!(prnt[1], cliq.index, TempBeliefMsg())
  end
  cda.upMsg = Dict{Symbol, BallTreeDensity}()
  cda.dwnMsg = Dict{Symbol, BallTreeDensity}()
  cda.upInitMsgs = Dict{Int, TempBeliefMsg}()
  cda.downInitMsg = TempBeliefMsg()
  setCliqStatus!(cliq, :null)
  setCliqDrawColor(cliq, "")
  return nothing
end

function resetCliqSolve!(dfg::G,
                         treel::BayesTree,
                         frt::Symbol;
                         solveKey::Symbol=:default  )::Nothing where G <: AbstractDFG
  #
  resetCliqSolve!(dfg, treel, getCliq(treel, frt), solveKey=solveKey)
end




"""
    $SIGNATURES

Inverse solve of predicted noise value and returns the associated "measured" noise value (also used as starting point for the solve).
"""
function solveFactorMeasurements(dfg::AbstractDFG,
                                 fctsym::Symbol  )
  #
  fcto = getFactor(dfg, fctsym)
  varsyms = fcto._variableOrderSymbols
  vars = map(x->getPoints(getKDE(dfg,x)), varsyms)
  fcttype = getFactorType(fcto)
  zDim = getData(fcto).fnc.zDim

  N = size(vars[1])[2]
  res = zeros(zDim)
  ud = FactorMetadata()
  meas = getSample(fcttype, N)
  meas0 = deepcopy(meas[1])

  function makemeas!(i, meas, dm)
    meas[1][:,i] = dm
    return meas
  end

  ggo = (i, dm) -> fcttype(res,ud,i,makemeas!(i, meas, dm),vars...)
  # ggo(1, [0.0;0.0])

  for idx in 1:N
    retry = 10
    while 0 < retry
      if isa(fcttype, FunctorPairwiseMinimize)
        r = optimize((x) -> ggo(idx,x), meas[1][:,idx]) # zeros(zDim)
        retry -= 1
        if !r.g_converged
          nsm = getSample(fcttype, 1)
          for count in 1:length(meas)
            meas[count][:,idx] = nsm[count][:,idx]
          end
        else
          break
        end
      elseif isa(fcttype, FunctorPairwise)
        ggnl = (rs, dm) -> fcttype(rs,ud,idx,makemeas!(idx, meas, dm),vars...)
        r = nlsolve(ggnl, meas[1][:,idx])
        break
      elseif isa(fcttype, FunctorSingleton)
        # assuming no partials at this point
        meas[1][:,:] .= vars[1][:,:]
        break
      end
    end
    # @assert meas[1][:,idx] == r.minimizer
  end

  # Gadfly.plot(z=(x,y)->ggo(1,[x;y]), xmin=[-pi],xmax=[pi],ymin=[-100.0],ymax=[100.0], Geom.contour)
  return meas[1], meas0
end


#
