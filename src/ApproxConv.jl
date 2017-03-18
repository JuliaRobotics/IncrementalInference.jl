
# currently monster of a spaghetti code mess, multiple types interacting at the same
# time. Safest at this point in development is to just get this code running and
# will refactor once the dust has settled.
# current code is the result of several unit tests between IIF and RoME.jl
# a unified test to follow -- after which refactoring can start
function evalPotentialSpecific{T <: FunctorPairwise}(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int64;
      N::Int64=100  )
  #
  # TODO -- enable partial constraints

  # TODO -- this part can be collapsed into common generic solver component, could be constructed and maintained at addFactor! time
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end
  # gwp.zDim[sfidx] ??
  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)
  # and return complete fr/gwp

  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  for gwp.particleidx in 1:maxlen
    # gwp(x, res)
    numericRootGenericRandomizedFnc!( fr )
        # r = nlsolve( gwp, ARR[sfidx][:,gwp.particleidx] )
        # remember this is a deepcopy of original sfidx, since we are generating a proposal distribution
        # and not directly replacing the existing variable belief estimate
        # gwp.params[gwp.varidx][:,gwp.particleidx] = r.zero[:]
  end

  return gwp.params[gwp.varidx]
  # return evalPotential(typ, Xi, solvefor)
end

function evalPotentialSpecific{T <: FunctorPairwiseMinimize}(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int64;
      N::Int64=100  )
  #
  # TODO -- this part can be collapsed into common generic solver component, could be constructed and maintained at addFactor! time
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end

  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)

  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  for gwp.particleidx in 1:maxlen
    # gwp(x, res)
    # implement minimization here
    res = zeros(fr.xDim)
    gg = (x) -> fr.gwp(x, res)
    r = optimize(  gg, fr.X[1:fr.xDim,fr.gwp.particleidx] )
    # TODO -- clearly lots of optmization to be done here
    fr.Y[1:fr.xDim] = r.minimizer
    fr.X[:,fr.gwp.particleidx] = fr.Y

    # error("not implemented yet")
  end
  return gwp.params[gwp.varidx]
end

function evalPotentialSpecific{T <: FunctorSingleton}(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      generalwrapper::GenericWrapParam{T},
      solvefor::Int64;
      N::Int64=100  )
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

function evalPotentialSpecific{T <: FunctorSingletonNH}(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      generalwrapper::GenericWrapParam{T},
      solvefor::Int64;
      N::Int64=100  )
  #

  # determine amount share of null hypothesis particles
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, N)

  # values of 0 imply null hypothesis
  # generalwrapper.usrfnc!.nullhypothesis::Distributions.Categorical
  nhc = rand(generalwrapper.usrfnc!.nullhypothesis, N) - 1

  val = getVal(Xi[1])
  d = size(val,1)
  var = Base.var(val,2)
  ENT = Distributions.MvNormal(zeros(d),10*diagm(var))

  for i in 1:N
    if nhc[i] == 0
      generalwrapper.measurement[:,i] = val[:,i] + rand(ENT)
    end
  end
  return generalwrapper.measurement
end

# function evalPotentialSpecific{T <: FunctorPartialSingleton}(
#       fnc::T,
#       Xi::Vector{Graphs.ExVertex},
#       generalwrapper::GenericWrapParam{T},
#       solvefor::Int64;
#       N::Int64=100  )
#   #
#   generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, N)
#   return generalwrapper.measurement[1]
# end

# Multiple dispatch occurs internally, resulting in factor graph potential evaluations
function evalFactor2(fgl::FactorGraph, fct::Graphs.ExVertex, solvefor::Int64; N::Int64=100)
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
    findRelatedFromPotential(fg, fctvert, varnodeid, N)

Compute proposal belief on varnodeid through fctvert representing some constraint in factor graph.
Always full dimension of variable node, where partial constraints will only influence directed
subset of variable dimensions. Remaining dimensions will keep existing variable values.
"""
function findRelatedFromPotential(fg::FactorGraph, idfct::Graphs.ExVertex, vertid::Int64, N::Int64) # vert
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
