
"""
    $(SIGNATURES)

Perform the nonlinear numerical operations to approximate the convolution with a particular user defined likelihood function (conditional), which as been prepared in the `frl` object.  This function uses root finding to enforce a non-linear function constraint.

Notes:
- remember this is a deepcopy of original sfidx, since we are generating a proposal distribution and not directly replacing the existing variable belief estimate

Future work:
- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
- improve handling of n and particleidx, especially considering future multithreading support

"""
function approxConvOnElements!(frl::FastRootGenericWrapParam{T},
                               elements::Union{Vector{Int}, UnitRange{Int}}  ) where {T <: FunctorPairwise}
  for n in elements
    frl.gwp.particleidx = n
    numericRootGenericRandomizedFnc!( frl )
  end
  # r = nlsolve( gwp, ARR[sfidx][:,gwp.particleidx] )
  # gwp.params[gwp.varidx][:,gwp.particleidx] = r.zero[:]
  nothing
end

"""
    $(SIGNATURES)

Perform the nonlinear numerical operations to approximate the convolution with a particular user defined likelihood function (conditional), which as been prepared in the `frl` object.  This function uses minimization of the res[1] variable.

Notes:
- remember this is a deepcopy of original sfidx, since we are generating a proposal distribution and not directly replacing the existing variable belief estimate
"""
function approxConvOnElements!(frl::FastRootGenericWrapParam{T},
                               elements::Union{Vector{Int}, UnitRange{Int}}) where {T <: FunctorPairwiseMinimize}
  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  res = zeros(frl.xDim)
  gg = (x) -> frl.gwp(res, x)
  for n in elements
    frl.gwp.particleidx = n
    res[:] = 0.0
    r = optimize( gg, frl.X[1:frl.xDim, frl.gwp.particleidx] )
    # TODO -- clearly lots of optmization to be done here
    frl.Y[1:frl.xDim] = r.minimizer
    frl.X[:,frl.gwp.particleidx] = frl.Y
  end
  nothing
end

"""
    $(SIGNATURES)

Prepare a common functor computation object `FastRootGenericWrapParam{T}` containing the `GenericWrapParam{T}` functor along with additional variables and information.
"""
function prepareFastRootGWP(gwp::GenericWrapParam{T},
                            Xi::Vector{Graphs.ExVertex},
                            solvefor::Int,
                            N::Int  ) where {T <: FunctorInferenceType}
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

"""
    $(SIGNATURES)

Common function to compute across a single user defined multi-hypothesis ambiguity per factor.  This function dispatches both `FunctorPairwise` and `FunctorPairwiseMinimize` factors.
"""
function computeAcrossHypothesis(frl::FastRootGenericWrapParam{T},
                                 allelements,
                                 activehypo,
                                 certainidx,
                                 sfidx) where {T <:Union{FunctorPairwise, FunctorPairwiseMinimize}}
  count = 0
  for (mhidx, vars) in activehypo
    count += 1
    if sfidx in certainidx || mhidx in certainidx # certainidx[count] in vars
      # standard case mhidx, sfidx = $mhidx, $sfidx
      frl.gwp.activehypo = vars
      approxConvOnElements!(frl, allelements[count])
    elseif mhidx == sfidx
      # multihypo, do conv case, mhidx == sfidx
      frl.gwp.activehypo = sort(union([sfidx;], certainidx))
      approxConvOnElements!(frl, allelements[count])
    elseif mhidx != sfidx
      # multihypo, take other value case
      # sfidx=2, mhidx=3:  2 should take a value from 3
      # sfidx=3, mhidx=2:  3 should take a value from 2
      frl.gwp.params[sfidx][:,allelements[count]] = view(frl.gwp.params[mhidx],:,allelements[count])
      # frl.gwp.params[sfidx][:,allelements[count]] = frl.gwp.params[mhidx][:,allelements[count]]
    else
      error("computeAcrossHypothesis -- not dealing with multi-hypothesis case correctly")
    end
  end
  nothing
end

"""
    $(SIGNATURES)

Prepare data required for null hypothesis cases during convolution.
"""
function assembleNullHypothesis(fr::FastRootGenericWrapParam{T},
                                maxlen::Int,
                                spreadfactor::Float64 ) where {T}
  #
  nhc = rand(fr.gwp.usrfnc!.nullhypothesis, maxlen) - 1
  val = fr.gwp.params[fr.gwp.varidx]
  d = size(val,1)
  var = Base.var(val,2) + 1e-3
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*diagm(var[:]))
  allelements = 1:maxlen
  return allelements, nhc, ENT
end

"""
    $(SIGNATURES)

Do true and null hypothesis computations based on data structures prepared earlier -- specific to `FunctorPairwiseNH`.  This function will be merged into a standard case for `FunctorPairwise/Minimize` in the future.
"""
function computeAcrossNullHypothesis!(frl::FastRootGenericWrapParam{T},
                                      allelements,
                                      nhc,
                                      ENT  ) where {T <: FunctorPairwiseNH}
  #
  # TODO --  Threads.@threads see area4 branch
  for n in allelements
    # frl.gwp(x, res)
    if nhc[n] != 0
      frl.gwp.particleidx = n
      numericRootGenericRandomizedFnc!( frl )
    else
      frl.gwp.params[frl.gwp.varidx][:,n] += rand(ENT)
    end
  end
  nothing
end


function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               gwp::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=100,
                               spreadfactor::Float64=10.0  ) where {T <: FunctorPairwiseNH}
  #

  # TODO -- could be constructed and maintained at addFactor! time
  fr, sfidx, maxlen = prepareFastRootGWP(gwp, Xi, solvefor, N)
  # prepare nullhypothesis
  allelements, nhc, ENT = assembleNullHypothesis(fr, maxlen, spreadfactor)

  # Compute across the true or null hypothesis
  computeAcrossNullHypothesis!(fr, allelements, nhc, ENT )

  return fr.gwp.params[gwp.varidx]
end

"""
    $(SIGNATURES)

Multiple dispatch wrapper for `<:FunctorPairwise` types, to prepare and execute the general approximate convolution with user defined factor residual functions.  This method also supports multihypothesis operations as one mechanism to introduce new modality into the proposal beliefs.
"""
function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               gwp::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=100  ) where {T <: Union{FunctorPairwise, FunctorPairwiseMinimize}}
  #
  fnc = gwp.usrfnc!

  # Prep computation variables
  fr, sfidx, maxlen = prepareFastRootGWP(gwp, Xi, solvefor, N)
  certainidx, allelements, activehypo, mhidx = assembleHypothesesElements!(fr.gwp.hypotheses, maxlen, sfidx, length(Xi))

  # perform the numeric solutions on the indicated elements
  computeAcrossHypothesis(fr, allelements, activehypo, certainidx, sfidx)

  return fr.gwp.params[gwp.varidx]
end


#  Singletons ==================================================================

"""
    $(SIGNATURES)

Multiple dispatch wrapper for evaluating the `genericwrapper::GenericWrapParam{<: FunctorSingleton}` types.
"""
function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               generalwrapper::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=0 ) where {T <: FunctorSingleton}
  #
  fnc = generalwrapper.usrfnc!

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

"""
    $(SIGNATURES)

Multiple dispatch wrapper for evaluating the `genericwrapper::GenericWrapParam{<: FunctorSingletonNH}` types.
Planned changes will fold null hypothesis in as a standard feature and no longer appear as a separate `InferenceType`.
"""
function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               generalwrapper::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=100,
                               spreadfactor::Float64=10.0 ) where {T <: FunctorSingletonNH}
  #
  fnc = generalwrapper.usrfnc!

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
      generalwrapper.measurement[1][:,i] = val[:,i] + rand(ENT)  # TODO use view and inplace add operation
    end
  end
  # TODO -- returning to memory location inside
  return generalwrapper.measurement[1]
end

"""
    $(SIGNATURES)

Single entry point for evaluating factors from factor graph, using multiple dispatch to locate the correct `evalPotentialSpecific` function.
"""
function evalFactor2(fgl::FactorGraph,
                     fct::Graphs.ExVertex,
                     solvefor::Int;
                     N::Int=100 )
  #

  gwp = getData(fct).fnc
  # TODO -- this build up of Xi is excessive and could happen at addFactor time
  Xi = Graphs.ExVertex[]
  count = 0
  for id in getData(fct).fncargvID
    count += 1
    xi = getVert(fgl,id)
    push!(Xi, xi ) # TODO localapi
    # push!(Xi, dlapi.getvertex(fgl,id))

    # TODO bad way to search for `solvefor`
    if xi.index == solvefor
      gwp.factormetadata.solvefor = Symbol(xi.label)
    end
  end
  return evalPotentialSpecific(Xi, gwp, solvefor, N=N)
end

"""
    $(SIGNATURES)

Draw samples from the approximate convolution of `towards` symbol using factor `fct` relative to the other variables.  In addition the `api` can be adjusted to recover the data from elsewhere (likely to be replaced/removed in the future).
"""
function approxConv(fgl::FactorGraph,
                    fct::Symbol,
                    towards::Symbol;
                    api::DataLayerAPI=localapi,
                    N::Int=-1  )
  #
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
