
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

```
# `allelements` example BearingRange [:x1, 0.5:l1a, 0.5:l1b]
# sfidx = (1=:x1,2=:l1a,3=:l1b)
if solvefor :x1, then allelem = [mhidx.==:l1a; mhidx.==l1b]
if solvefor :l1a, then allelem = [mhidx.==:l1a] and ARR[solvefor][:,mhidx.==:l1b]=ARR[:l1b][:,mhidx.==:l1b]
if solvefor :l1b, then allelem = [mhidx.==:l1b] and ARR[solvefor][:,mhidx.==:l1a]=ARR[:l1a][:,mhidx.==:l1a]
if solvefor 1, then allelem = [mhidx.==2; mhidx.==3]
if solvefor 2, then allelem = [mhidx.==2] and ARR[solvefor][:,mhidx.==3]=ARR[3][:,mhidx.==3]
if solvefor 3, then allelem = [mhidx.==3] and ARR[solvefor][:,mhidx.==2]=ARR[2][:,mhidx.==2]

# `activehypo` in example mh=[0;0.5;0.5]
sfidx=1, mhidx=2:  ah = [1;2]
sfidx=1, mhidx=3:  ah = [1;3]
sfidx=2, mhidx=2:  ah = [1;2]
sfidx=2, mhidx=3:  ah = # 2 should take a value from 3
sfidx=3, mhidx=2:  ah = # 3 should take a value from 2
sfidx=3, mhidx=3:  ah = [1;3]

# `activehypo` in example mh=[0;0.33;0.33;0.34]
sfidx=1, mhidx=2:  ah = [1;2]
sfidx=1, mhidx=3:  ah = [1;3]
sfidx=1, mhidx=4:  ah = [1;4]
...
```
"""
function assembleHypothesesElements!(allelements::Array,
            activehypo::Array,
            mh::Categorical,
            maxlen::Int,
            sfidx::Int,
            mhidx,
            lenXi::Int  )
  #
  # @show mhidx
  allidx = 1:maxlen

  # this is not going to work? sfidx could be anything
  if mh.p[sfidx] < 1e-10
    pidx = 0
    for pval in mh.p
      pidx += 1
      if pval > 1e-10
        iterarr = allidx[mhidx .== pidx]
        push!(allelements, iterarr)
        iterah = sort([sfidx;pidx]) # TODO -- currently only support binary factors in multihypo mode
        push!(activehypo, iterah)
      end
    end
  elseif mh.p[sfidx] > 1e-10
    pidx = 0
    for pval in mh.p
      pidx += 1
      if pval > 1e-10 && sfidx == pidx
        iterarr = allidx[mhidx .== pidx]
        push!(allelements, iterarr)
        iterah = sort([sfidx;pidx]) # TODO -- currently only support binary factors in multihypo mode
        push!(activehypo, iterah)
      end
    end
    # must still include cases where sfidx != pidx
  else
    error("Unknown hypothesis case, got sfidx=$(sfidx) with mh.p=$(mh.p)")
  end

  nothing
end
function assembleHypothesesElements!(allelements::Array, activehypo::Array, mh::Void, maxlen::Int, sfidx::Int, mhidx, lenXi::Int)
  # error("assembleHypothesesElements!(..) -- Error in code design, refactor of general multihypothesis situations required if you arrived here.")
  push!(allelements, 1:maxlen)
  push!(activehypo, 1:lenXi)
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
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end

  # TODO -- introduce special case for multihypothesis
  assembleHypothesesElements!(allelements, activehypo, gwp.hypotheses, maxlen, sfidx, mhidx, length(Xi))
  @show size(allelements), size(activehypo)
  # if gwp.hypotheses == nothing
  #   push!(allelements, 1:maxlen)
  # else
  #   # consider merging with prepareparamsarray
  # end


  # Construct complete fr (with fr.gwp) object
  # TODO -- create FastRootGenericWrapParam at addFactor time only
  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)

  # gwp.activehypo = 1:length(Xi)# for regular n-ary operation
  # different hypotheses require

   ##if have .p# mhwhoszero = findfirst(gwp.hypotheses.p .< 1e-10)
  # perform the numeric solutions on the indicated elements
  for idx in 1:length(allelements)
    @show gwp.activehypo = activehypo[idx]
    approxConvOnElements!(fr, allelements[idx])
  end
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
