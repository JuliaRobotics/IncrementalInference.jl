"""
    $(SIGNATURES)

Perform the nonlinear numerical operations to approximate the convolution with a particular user defined likelihood function (conditional), which as been prepared in the `frl` object.  This function uses root finding to enforce a non-linear function constraint.

Notes:
- remember this is a deepcopy of original sfidx, since we are generating a proposal distribution and not directly replacing the existing variable belief estimate

Future work:
- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
- improve handling of n and particleidx, especially considering future multithreading support

"""
function approxConvOnElements!(ccwl::CommonConvWrapper{T},
                               elements::Union{Vector{Int}, UnitRange{Int}}, ::Type{MultiThreaded}  ) where {T <: Union{FunctorPairwise, FunctorPairwiseMinimize}}
  Threads.@threads for n in elements
    # ccwl.thrid_ = Threads.threadid()
    ccwl.cpt[Threads.threadid()].particleidx = n
    # ccall(:jl_, Nothing, (Any,), "starting loop, thrid_=$(Threads.threadid()), partidx=$(ccwl.cpt[Threads.threadid()].particleidx)")
    numericRootGenericRandomizedFnc!( ccwl )
  end
  nothing
end

"""
    $(SIGNATURES)

Perform the nonlinear numerical operations to approximate the convolution with a particular user defined likelihood function (conditional), which as been prepared in the `frl` object.  This function uses root finding to enforce a non-linear function constraint.

Notes:
- remember this is a deepcopy of original sfidx, since we are generating a proposal distribution and not directly replacing the existing variable belief estimate

Future work:
- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
- improve handling of n and particleidx, especially considering future multithreading support

"""
function approxConvOnElements!(ccwl::CommonConvWrapper{T},
                               elements::Union{Vector{Int}, UnitRange{Int}}, ::Type{SingleThreaded}) where {T <: Union{FunctorPairwise, FunctorPairwiseMinimize}}
  #
  for n in elements
    ccwl.cpt[Threads.threadid()].particleidx = n
    numericRootGenericRandomizedFnc!( ccwl )
  end
  nothing
end


"""
    $(SIGNATURES)

Perform the nonlinear numerical operations to approximate the convolution with a particular user defined likelihood function (conditional), which as been prepared in the `frl` object.  This function uses root finding to enforce a non-linear function constraint.

Notes:
- remember this is a deepcopy of original sfidx, since we are generating a proposal distribution and not directly replacing the existing variable belief estimate

Future work:
- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
- improve handling of n and particleidx, especially considering future multithreading support

"""
function approxConvOnElements!(ccwl::CommonConvWrapper{T},
                               elements::Union{Vector{Int}, UnitRange{Int}} )  where {T <: Union{FunctorPairwise, FunctorPairwiseMinimize}}
  #
  approxConvOnElements!(ccwl, elements, ccwl.threadmodel)
end



"""
    $(SIGNATURES)

Prepare a common functor computation object `prepareCommonConvWrapper{T}` containing the user factor functor along with additional variables and information using during approximate convolution computations.
"""
function prepareCommonConvWrapper!(ccwl::CommonConvWrapper{T},
                                   Xi::Vector{DFGVariable},
                                   solvefor::Symbol,
                                   N::Int  ) where {T <: FunctorInferenceType}
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
  # should be selecting for the correct multihypothesis mode here with `gwp.params=ARR[??]`
  ccwl.params = ARR
  ccwl.measurement = getSample(ccwl.usrfnc!, maxlen) # ccwl.samplerfnc
  if ccwl.specialzDim
    ccwl.zDim = ccwl.usrfnc!.zDim[sfidx]
  else
    ccwl.zDim = size(ccwl.measurement[1],1) # TODO -- zDim aspect needs to be reviewed
  end
  ccwl.varidx = sfidx

  ccwl.xDim = size(ccwl.params[sfidx],1)
  # ccwl.xDim = size(ccwl.cpt[1].X,1)
  # info("what? sfidx=$(sfidx), ccwl.xDim = size(ccwl.params[sfidx]) = $(ccwl.xDim), size=$(size(ccwl.params[sfidx]))")
  for i in 1:Threads.nthreads()
    ccwl.cpt[i].X = ccwl.params[sfidx]
    ccwl.cpt[i].p = collect(1:size(ccwl.cpt[i].X,1))
    ccwl.cpt[i].Y = zeros(ccwl.xDim)
    ccwl.cpt[i].res = zeros(ccwl.xDim) # used in ccw functor for FunctorPairwiseMinimize
  end

  return sfidx, maxlen
end

"""
    $(SIGNATURES)

Common function to compute across a single user defined multi-hypothesis ambiguity per factor.  This function dispatches both `FunctorPairwise` and `FunctorPairwiseMinimize` factors.
"""
function computeAcrossHypothesis!(ccwl::CommonConvWrapper{T},
                                  allelements,
                                  activehypo,
                                  certainidx,
                                  sfidx) where {T <:Union{FunctorPairwise, FunctorPairwiseMinimize}}
  count = 0
  # TODO remove assert once all GenericWrapParam has been removed
  # @assert norm(ccwl.certainhypo - certainidx) < 1e-6
  for (mhidx, vars) in activehypo
    count += 1
    if sfidx in certainidx || mhidx in certainidx || mhidx == sfidx
      # hypo case mhidx, sfidx = $mhidx, $sfidx
      for i in 1:Threads.nthreads()  ccwl.cpt[i].activehypo = vars; end
      approxConvOnElements!(ccwl, allelements[count])
    # elseif mhidx == sfidx
    #   # multihypo, do conv case, mhidx == sfidx
    #   ah = sort(union([sfidx;], certainidx))
    #   @assert norm(ah - vars) < 1e-10
    #   for i in 1:Threads.nthreads()  ccwl.cpt[i].activehypo = ah; end
    #   approxConvOnElements!(ccwl, allelements[count])
    elseif mhidx != sfidx   # sfidx in uncertnidx
      # multihypo, take other value case
      # sfidx=2, mhidx=3:  2 should take a value from 3
      # sfidx=3, mhidx=2:  3 should take a value from 2
      ccwl.params[sfidx][:,allelements[count]] = view(ccwl.params[mhidx],:,allelements[count])
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
function assembleNullHypothesis(ccwl::CommonConvWrapper{T},
                                maxlen::Int,
                                spreadfactor::Float64 ) where {T}
  #
  nhc = rand(ccwl.usrfnc!.nullhypothesis, maxlen) .- 1
  val = ccwl.params[ccwl.varidx]
  d = size(val,1)
  var = Statistics.var(val,dims=2) .+ 1e-3
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*Matrix(Diagonal(var[:])))
  allelements = 1:maxlen
  return allelements, nhc, ENT
end

"""
    $(SIGNATURES)

Do true and null hypothesis computations based on data structures prepared earlier -- specific to `FunctorPairwiseNH`.  This function will be merged into a standard case for `FunctorPairwise/Minimize` in the future.
"""
function computeAcrossNullHypothesis!(ccwl::CommonConvWrapper{T},
                                      allelements,
                                      nhc,
                                      ENT  ) where {T <: FunctorPairwiseNH}
  #
  # TODO --  Threads.@threads see area4 branch
  for n in allelements
    # ccwl.gwp(x, res)
    if nhc[n] != 0
      ccwl.cpt[Threads.threadid()].particleidx = n
      numericRootGenericRandomizedFnc!( ccwl )
    else
      ccwl.params[ccwl.varidx][:,n] += rand(ENT)
    end
  end
  nothing
end



"""
    $(SIGNATURES)

Multiple dispatch wrapper for `<:FunctorPairwise` types, to prepare and execute the general approximate convolution with user defined factor residual functions.  This method also supports multihypothesis operations as one mechanism to introduce new modality into the proposal beliefs.

Planned changes will fold null hypothesis in as a standard feature and no longer appear as a separate `InferenceType`.
"""
function evalPotentialSpecific(Xi::Vector{DFGVariable},
                               ccwl::CommonConvWrapper{T},
                               solvefor::Symbol;
                               N::Int=100,
                               spreadfactor::Float64=10.0,
                               dbg::Bool=false ) where {T <: FunctorPairwiseNH}
  #

  # TODO -- could be constructed and maintained at addFactor! time
  sfidx, maxlen = prepareCommonConvWrapper!(ccwl, Xi, solvefor, N)
  # prepare nullhypothesis
  allelements, nhc, ENT = assembleNullHypothesis(ccwl, maxlen, spreadfactor)

  # Compute across the true or null hypothesis
  computeAcrossNullHypothesis!(ccwl, allelements, nhc, ENT )

  return ccwl.params[ccwl.varidx]
end

function evalPotentialSpecific(Xi::Vector{DFGVariable},
                               ccwl::CommonConvWrapper{T},
                               solvefor::Symbol;
                               N::Int=100,
                               dbg::Bool=false ) where {T <: Union{FunctorPairwise, FunctorPairwiseMinimize}}
  #

  # Prep computation variables
  sfidx, maxlen = prepareCommonConvWrapper!(ccwl, Xi, solvefor, N)
  _, allelements, activehypo, mhidx = assembleHypothesesElements!(ccwl.hypotheses, maxlen, sfidx, length(Xi))
  certainidx = ccwl.certainhypo

  # perform the numeric solutions on the indicated elements
  # error("ccwl.xDim=$(ccwl.xDim)")
  computeAcrossHypothesis!(ccwl, allelements, activehypo, certainidx, sfidx)

  return ccwl.params[ccwl.varidx]
end

function evalPotentialSpecific(Xi::Vector{DFGVariable},
                               ccwl::CommonConvWrapper{T},
                               solvefor::Symbol;
                               N::Int=0,
                               dbg::Bool=false ) where {T <: FunctorSingleton}
  #
  fnc = ccwl.usrfnc!

  nn = (N <= 0 ? size(getVal(Xi[1]),2) : N)
  ccwl.measurement = getSample(ccwl.usrfnc!, nn)
  if !ccwl.partial
    return ccwl.measurement[1]
  else
    val = deepcopy(getVal(Xi[1]))
    i = 0
    for dimnum in fnc.partial
      i += 1
      val[dimnum,:] = ccwl.measurement[1][i,:]
    end
    return val
  end
end

function evalPotentialSpecific(Xi::Vector{DFGVariable},
                               ccwl::CommonConvWrapper{T},
                               solvefor::Symbol;
                               N::Int=100,
                               spreadfactor::Float64=10.0,
                               dbg::Bool=false ) where {T <: FunctorSingletonNH}
  #
  fnc = ccwl.usrfnc!

  val = getVal(Xi[1])
  d = size(val,1)
  var = Statistics.var(val, dims=2) .+ 1e-3

  # determine amount share of null hypothesis particles
  ccwl.measurement = getSample(ccwl.usrfnc!, N)
  # values of 0 imply null hypothesis
  # ccwl.usrfnc!.nullhypothesis::Distributions.Categorical
  nhc = rand(ccwl.usrfnc!.nullhypothesis, N) .- 1

  # TODO -- not valid for manifold
  # TODO bad memory management
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*Matrix(Diagonal(var[:])) )

  for i in 1:N
    if nhc[i] == 0
      ccwl.measurement[1][:,i] = val[:,i] + rand(ENT)  # TODO use view and inplace add operation
    end
  end
  # TODO -- returning to memory location inside
  return ccwl.measurement[1]
end

"""
    $(SIGNATURES)

Single entry point for evaluating factors from factor graph, using multiple dispatch to locate the correct `evalPotentialSpecific` function.
"""
function evalFactor2(dfg::G,
                     fct::DFGFactor,
                     solvefor::Symbol;
                     N::Int=100,
                     dbg::Bool=false) where G <: AbstractDFG
  #

  ccw = solverData(fct).fnc
  # TODO -- this build up of Xi is excessive and could happen at addFactor time
  Xi = DFGVariable[]
  count = 0
  # TODO replace fncargvID with neighbors
  variablelist = Vector{Symbol}(undef, length(solverData(fct).fncargvID))
  for id in solverData(fct).fncargvID
    count += 1
    xi = DFG.getVariable(dfg, id)
    push!(Xi, xi )
    # push!(Xi, dlapi.getvertex(fgl,id))

    # TODO do only once at construction time -- staring it here to be sure the code is calling factors correctly
    variablelist[count] = xi.label

    # TODO bad way to search for `solvefor`
    if xi.label == solvefor
      for i in 1:Threads.nthreads()
        ccw.cpt[i].factormetadata.solvefor = xi.label
      end
    end
  end
  for i in 1:Threads.nthreads()
    ccw.cpt[i].factormetadata.variablelist = variablelist
  end
  return evalPotentialSpecific(Xi, ccw, solvefor, N=N, dbg=dbg)
end
# import IncrementalInference: evalFactor2, approxConv
"""
    $(SIGNATURES)

Draw samples from the approximate convolution of `towards` symbol using factor `fct` relative to the other variables.  In addition the `api` can be adjusted to recover the data from elsewhere (likely to be replaced/removed in the future).
"""
function approxConv(dfg::G,
                    fct::Symbol,
                    towards::Symbol;
                    N::Int=-1  ) where G <: AbstractDFG
  #
  fc = getFactor(dfg, fct)
  v1 = getVariable(dfg, towards)
  N = N == -1 ? getNumPts(v1) : N
  return evalFactor2(dfg, fc, v1.label, N=N)
end


function approxConvBinary(arr::Array{Float64,2}, meas::T, outdims::Int; N::Int=0) where {T <: FunctorInferenceType}
  N = N == 0 ? size(arr,2) : N
  pts = zeros(outdims,N);
  t = Array{Array{Float64,2},1}()
  push!(t,arr)
  push!(t,pts)

  measurement = getSample(meas, N)

  ccw = CommonConvWrapper(meas, t[2], size(measurement[1],2), t, varidx=2, measurement=measurement)

  for n in 1:N
    ccw.cpt[Threads.threadid()].particleidx = n
    numericRootGenericRandomizedFnc!( ccw )
  end
  return pts
end



"""
    $(SIGNATURES)

Compute proposal belief on `vertid` through `fct` representing some constraint in factor graph.
Always full dimension variable node -- partial constraints will only influence subset of variable dimensions.
The remaining dimensions will keep pre-existing variable values.

Notes
- fulldim is true when "rank-deficient" -- TODO swap to false (or even float)
"""
function findRelatedFromPotential(dfg::G,
                                  fct::DFGFactor,
                                  varid::Symbol,
                                  N::Int,
                                  dbg::Bool=false  )::Tuple{BallTreeDensity,Float64} where G <: AbstractDFG
  # assuming it is properly initialized TODO
  ptsbw = evalFactor2(dfg, fct, varid, N=N, dbg=dbg);
  # determine if evaluation is "dimension-deficient"

  # solvable dimension
  inferdim = getFactorSolvableDim(dfg, fct, varid)
  # zdim = getFactorDim(fct)
  # vdim = getVariableDim(DFG.getVariable(dfg, varid))

  # TODO -- better to upsample before the projection
  Ndim = size(ptsbw,1)
  Npoints = size(ptsbw,2)
  # Assume we only have large particle population sizes, thanks to addNode!
  manis = getManifolds(dfg, varid)
  # manis = getSofttype(DFG.getVariable(dfg, varid)).manifolds # older
  p = AMP.manikde!(ptsbw, manis)
  if Npoints != N # this is where we control the overall particle set size
      p = resample(p,N)
  end
  return p, inferdim
end




#
