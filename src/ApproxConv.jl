
"""
    $(SIGNATURES)

Perform the nonlinear numerical operations to approximate the convolution with a particular user defined likelihood function (conditional), which as been prepared in the `frl` object.
"""
function approxConvOnElements!(frl::FastRootGenericWrapParam, elements::Union{Vector{Int}, UnitRange{Int}})
  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  # TODO -- improve handling of n and particleidx, especially considering future multithreading support
  for n in elements
    frl.gwp.particleidx = n
    numericRootGenericRandomizedFnc!( frl )
  end
  # r = nlsolve( gwp, ARR[sfidx][:,gwp.particleidx] )
  # remember this is a deepcopy of original sfidx, since we are generating a proposal distribution
  # and not directly replacing the existing variable belief estimate
  # gwp.params[gwp.varidx][:,gwp.particleidx] = r.zero[:]
  nothing
end

function prepareFastRootGWP(T::Type, gwp, Xi::Vector{Graphs.ExVertex}, solvefor::Int, N::Int)
  # TODO -- this part can be collapsed into common generic solver component
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
  # should be selecting for the correct multihypothesis mode here with `gwp.params=ARR[??]`
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  size(gwp.measurement[1])
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end
  # Construct complete fr (with fr.gwp) object
  # TODO -- create FastRootGenericWrapParam at addFactor time only?
  FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp), sfidx, maxlen
end

function computeAcrossHypothesis(fr, allelements, activehypo, certainidx, sfidx)
  count = 0
  for (mhidx, vars) in activehypo
    count += 1
    if sfidx in certainidx || mhidx in certainidx # certainidx[count] in vars
      # standard case mhidx, sfidx = $mhidx, $sfidx
      fr.gwp.activehypo = vars
      approxConvOnElements!(fr, allelements[count])
    elseif mhidx == sfidx
      # multihypo, do conv case, mhidx == sfidx
      fr.gwp.activehypo = sort(union([sfidx;], certainidx))
      approxConvOnElements!(fr, allelements[count])
    elseif mhidx != sfidx
      # multihypo, take other value case
      # sfidx=2, mhidx=3:  2 should take a value from 3
      # sfidx=3, mhidx=2:  3 should take a value from 2
      fr.gwp.params[sfidx][:,allelements[count]] = fr.gwp.params[mhidx][:,allelements[count]]
    else
      error("evalPotentialSpecific -- not dealing with multi-hypothesis case correctly")
    end
  end
  nothing
end

"""
    $(SIGNATURES)

Multiple dispatch wrapper for `<:FunctorPairwise` types, to prepare and execute the general approximate convolution with user defined factor residual functions.  This method also supports multihypothesis operations as one mechanism to introduce new modality into the proposal beliefs.
"""
function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int;
      N::Int=100  ) where {T <: FunctorPairwise}
  #

  # Prep computation variables
  fr, sfidx, maxlen = prepareFastRootGWP(T, gwp, Xi, solvefor, N)
  certainidx, allelements, activehypo, mhidx = assembleHypothesesElements!(fr.gwp.hypotheses, maxlen, sfidx, length(Xi))

  # perform the numeric solutions on the indicated elements
  computeAcrossHypothesis(fr, allelements, activehypo, certainidx, sfidx)

  return fr.gwp.params[gwp.varidx]
end


function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int;
      N::Int=100,
      spreadfactor::Float64=10.0  ) where {T <: FunctorPairwiseNH}
  #

  # TODO -- could be constructed and maintained at addFactor! time
  fr, sfidx, maxlen = prepareFastRootGWP(T, gwp, Xi, solvefor, N)

  # nullhypothesis
  nhc = rand(fr.gwp.usrfnc!.nullhypothesis, maxlen) - 1
  val = fr.gwp.params[fr.gwp.varidx]
  d = size(val,1)
  var = Base.var(val,2) + 1e-3
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*diagm(var[:]))

  allelements = 1:maxlen

  # TODO --  Threads.@threads see area4 branch
  for n in allelements
    # fr.gwp(x, res)
    if nhc[n] != 0
      fr.gwp.particleidx = n
      numericRootGenericRandomizedFnc!( fr )
    else
      fr.gwp.params[fr.gwp.varidx][:,n] += rand(ENT)
    end
  end

  return fr.gwp.params[gwp.varidx]
end


function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int;
      N::Int=100  ) where {T <: FunctorPairwiseMinimize}
  #
  # TODO -- could be constructed and maintained at addFactor! time
  fr, sfidx, maxlen = prepareFastRootGWP(T, gwp, Xi, solvefor, N)

  allelements = 1:maxlen

  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  for n in allelements
    fr.gwp.particleidx = n
    res = zeros(fr.xDim)
    gg = (x) -> fr.gwp(res, x)
    r = optimize( gg, fr.X[1:fr.xDim, fr.gwp.particleidx] )
    # TODO -- clearly lots of optmization to be done here
    fr.Y[1:fr.xDim] = r.minimizer
    fr.X[:,fr.gwp.particleidx] = fr.Y
  end
  return fr.gwp.params[gwp.varidx]
end


#  Singletons ==================================================================

function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      generalwrapper::GenericWrapParam{T},
      solvefor::Int;
      N::Int=0  ) where {T <: FunctorSingleton}
  #
  nn = N != 0 ? N : size(getVal(Xi[1]),2)
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, nn)
  if !generalwrapper.partial
    return generalwrapper.measurement[1]
  else
    val = deepcopy(getVal(Xi[1]))
    i = 0
    for dimnum in fnc.partial
      i += 1
      val[dimnum,:] = generalwrapper.measurement[1][i,:]
    end
    return val
  end
end

function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      generalwrapper::GenericWrapParam{T},
      solvefor::Int;
      N::Int=100,
      spreadfactor::Float64=10.0  ) where {T <: FunctorSingletonNH}
  #
  val = getVal(Xi[1])
  d = size(val,1)
  var = Base.var(val,2) + 1e-3

  # determine amount share of null hypothesis particles
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, N)
  # values of 0 imply null hypothesis
  # generalwrapper.usrfnc!.nullhypothesis::Distributions.Categorical
  nhc = rand(generalwrapper.usrfnc!.nullhypothesis, N) - 1

  # TODO -- not valid for manifold
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*diagm(var[:]))

  for i in 1:N
    if nhc[i] == 0
      generalwrapper.measurement[1][:,i] = val[:,i] + rand(ENT)
    end
  end
  # TODO -- returning to memory location inside
  return generalwrapper.measurement[1]
end


# Multiple dispatch occurs internally, resulting in factor graph potential evaluations
function evalFactor2(fgl::FactorGraph, fct::Graphs.ExVertex, solvefor::Int; N::Int=100)
  # return evalPotential(fct.attributes["data"].fnc, solvefor) #evalPotential(fct.attributes["fnc"], solvefor)

  # TODO -- this build up of Xi is excessive and should be reduced
  # could happen at addFactor time
  Xi = Graphs.ExVertex[]
  for id in getData(fct).fncargvID
    push!(Xi, getVert(fgl,id) ) # TODO localapi
    # push!(Xi, dlapi.getvertex(fgl,id))
  end
  fnctype = getData(fct).fnc
  # fnctype = fct.attributes["data"].fnc
  return evalPotentialSpecific(fnctype.usrfnc!, Xi, fnctype, solvefor, N=N)
  # return evalPotentialSpecific(modulefnc, Xi, fnctype, solvefor, N=N)
end

function approxConv(fgl::FactorGraph, fct::Symbol, towards::Symbol; api::DataLayerAPI=localapi, N=-1)
  fc = getVert(fgl, fct, nt=:fct, api=api)
  v1 = getVert(fgl, towards, api=api)
  N = N == -1 ? N : getNumPts(v1)
  return evalFactor2(fgl, fc, v1.index, N=N)
end


"""
    $(SIGNATURES)

Compute proposal belief on varnodeid through fctvert representing some constraint in factor graph.
Always full dimension of variable node, where partial constraints will only influence directed
subset of variable dimensions. Remaining dimensions will keep existing variable values.
"""
function findRelatedFromPotential(fg::FactorGraph, idfct::Graphs.ExVertex, vertid::Int, N::Int) # vert
  # assuming it is properly initialized TODO
  ptsbw = evalFactor2(fg, idfct, vertid, N=N);
  # sum(abs(ptsbw)) < 1e-14 ? error("findRelatedFromPotential -- an input is zero") : nothing  # NOTE -- disable this validation test

  # TODO -- better to upsample before the projection
  Ndim = size(ptsbw,1)
  Npoints = size(ptsbw,2)
  # Assume we only have large particle population sizes, thanks to addNode!
  p = kde!(ptsbw, "lcv")
  if Npoints != N # this is where we control the overall particle set size
      p = resample(p,N)
  end
  return p
end
