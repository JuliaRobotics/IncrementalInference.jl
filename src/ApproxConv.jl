
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
  # currently monster of a spaghetti code mess (WIP) with multiple types interacting at the same
  # time.  Safest is to ensure code is producing correct results and refactor with unit tests in place.

  allelements = []
  activehypo = []
  # TODO -- this part can be collapsed into common generic solver component
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx, mhidx = prepareparamsarray!(ARR, Xi, N, solvefor, gwp.hypotheses)
  # should be selecting for the correct multihypothesis mode here with `gwp.params=ARR[??]`
  gwp.params = ARR
  @show gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end

  # TODO -- introduce special case for multihypothesis
  certainidx = assembleHypothesesElements!(allelements, activehypo, gwp.hypotheses, maxlen, sfidx, mhidx, length(Xi))
  # @show size(allelements), size(activehypo), activehypo

  # Construct complete fr (with fr.gwp) object
  # TODO -- create FastRootGenericWrapParam at addFactor time only
  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)

  # perform the numeric solutions on the indicated elements
  # @show sfidx
  # @show certainidx
  # @show activehypo
  @show allelements
  # mhiters = gwp.hypotheses == nothing ? length(allelements) : length(gwp.hypotheses.p)
  # normalidx = 1

  count = 0
  for (mhidx, vars) in activehypo
    # @show mhidx, vars
    @show count += 1
    @show sfidx, mhidx, vars, certainidx, count
    @show length(allelements[count])
    if sfidx in certainidx || mhidx in certainidx # certainidx[count] in vars
      @show gwp.activehypo = vars
      approxConvOnElements!(fr, allelements[count])
    elseif mhidx == sfidx
      # nothing to be done
      @show mvars = sort(union([sfidx;], certainidx))
      @show gwp.activehypo = mvars
      approxConvOnElements!(fr, allelements[count])
      info("multihypo, do conv case")
    # elseif length(allelements[count]) > 0
    #   @show gwp.activehypo = vars
    #   approxConvOnElements!(fr, allelements[count])
    else
      error("deal with mh case")
    end
  end

  # TODO -- refactor to always use n-number of Discrete Categoricals depending on how many the user wants to incorporate
  # currently marginalizing only one generalized discrete variable under multihypo interface (simplification for agnostic user)
  # for idx in 1:mhiters
  #   @show idx, activehypo, activehypo[normalidx][2]
  #   @show certainidx[1] in activehypo[normalidx][2]
  #   if idx == activehypo[normalidx][1] && certainidx[1] in activehypo[normalidx][2]
  #     # general case
  #     #!(idx in certainidx) && certainidx[1] in activehypo[normalidx][2]
  #     @show normalidx
  #     @show gwp.activehypo = activehypo[normalidx][2]
  #     approxConvOnElements!(fr, allelements[normalidx])
  #     @show normalidx += 1
  #   else
  #     # deal with multihypothesis cases
  #     @show sfidx, idx, activehypo[idx], size(mhidx)
  #     error("deal with mh case")
  #   # else
  #   #   error("Unknown multihypothesis case")
  #   end
  # end
  # if !(mhwhoszero in activehypo[idx])

  return gwp.params[gwp.varidx]
end


function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int;
      N::Int=100,
      spreadfactor::Float64=10.0  ) where {T <: FunctorPairwiseNH}
  #

  # TODO -- this part can be collapsed into common generic solver component, could be constructed and maintained at addFactor! time
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx, mhidx = prepareparamsarray!(ARR, Xi, N, solvefor, gwp.hypotheses)
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end
  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)
  # and return complete fr/gwp

  # nullhypothesis
  nhc = rand(gwp.usrfnc!.nullhypothesis, maxlen) - 1
  val = gwp.params[gwp.varidx]
  d = size(val,1)
  var = Base.var(val,2) + 1e-3
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*diagm(var[:]))

  allelements = 1:maxlen

  # TODO --  Threads.@threads see area4 branch
  for n in allelements
    # gwp(x, res)
    if nhc[n] != 0
      gwp.particleidx = n
      numericRootGenericRandomizedFnc!( fr )
    else
      gwp.params[gwp.varidx][:,n] += rand(ENT)
    end
  end

  return gwp.params[gwp.varidx]
  # return evalPotential(typ, Xi, solvefor)
end


function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int;
      N::Int=100  ) where {T <: FunctorPairwiseMinimize}
  #
  # TODO -- this part can be collapsed into common generic solver component, could be constructed and maintained at addFactor! time
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx, mhidx = prepareparamsarray!(ARR, Xi, N, solvefor, gwp.hypotheses)
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end

  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)

  allelements = 1:maxlen

  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  for n in allelements
    gwp.particleidx = n
    # gwp(x, res)
    # implement minimization here
    res = zeros(fr.xDim)
    gg = (x) -> fr.gwp(res, x)
    r = optimize(  gg, fr.X[1:fr.xDim,fr.gwp.particleidx] )
    # TODO -- clearly lots of optmization to be done here
    fr.Y[1:fr.xDim] = r.minimizer
    fr.X[:,fr.gwp.particleidx] = fr.Y

    # error("not implemented yet")
  end
  return gwp.params[gwp.varidx]
end

function evalPotentialSpecific(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      generalwrapper::GenericWrapParam{T},
      solvefor::Int;
      N::Int=100  ) where {T <: FunctorSingleton}
  #
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, N)
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

  # determine amount share of null hypothesis particles
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, N)
  # values of 0 imply null hypothesis
  # generalwrapper.usrfnc!.nullhypothesis::Distributions.Categorical
  nhc = rand(generalwrapper.usrfnc!.nullhypothesis, N) - 1

  val = getVal(Xi[1])
  d = size(val,1)
  var = Base.var(val,2) + 1e-3
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
  # for id in fct.attributes["data"].fncargvID
    push!(Xi, getVert(fgl,id) ) # TODO localapi
    # push!(Xi, dlapi.getvertex(fgl,id))
  end
  fnctype = getData(fct).fnc
  # fnctype = fct.attributes["data"].fnc
  return evalPotentialSpecific(fnctype.usrfnc!, Xi, fnctype, solvefor, N=N)
  # return evalPotentialSpecific(modulefnc, Xi, fnctype, solvefor, N=N)
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
