
export calcFactorResidual
export findRelatedFromPotential



"""
    $(SIGNATURES)

Perform the nonlinear numerical operations to approximate the convolution with a particular user defined likelihood function (conditional), which as been prepared in the `frl` object.  This function uses root finding to enforce a non-linear function constraint.

Notes:
- remember this is a deepcopy of original sfidx, since we are generating a proposal distribution and not directly replacing the existing variable belief estimate

Future work:
- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
- improve handling of n and particleidx, especially considering future multithreading support
"""
function approxConvOnElements!( ccwl::Union{CommonConvWrapper{F},
                                            CommonConvWrapper{Mixture{N_,F,S,T}}},
                                elements::Union{Vector{Int}, UnitRange{Int}}, ::Type{MultiThreaded}  ) where {N_,F<:AbstractRelative,S,T}
  #
  Threads.@threads for n in elements
    # ccwl.thrid_ = Threads.threadid()
    ccwl.cpt[Threads.threadid()].particleidx = n
    
    # ccall(:jl_, Nothing, (Any,), "starting loop, thrid_=$(Threads.threadid()), partidx=$(ccwl.cpt[Threads.threadid()].particleidx)")
    _solveCCWNumeric!( ccwl )
  end
  nothing
end


function approxConvOnElements!( ccwl::Union{CommonConvWrapper{F},
                                            CommonConvWrapper{Mixture{N_,F,S,T}}},
                                elements::Union{Vector{Int}, UnitRange{Int}}, ::Type{SingleThreaded}) where {N_,F<:AbstractRelative,S,T}
  #
  for n in elements
    ccwl.cpt[Threads.threadid()].particleidx = n    
    _solveCCWNumeric!( ccwl )
  end
  nothing
end


function approxConvOnElements!( ccwl::Union{CommonConvWrapper{F},
                                            CommonConvWrapper{Mixture{N_,F,S,T}}},
                                elements::Union{Vector{Int}, UnitRange{Int}} )  where {N_,F<:AbstractRelative,S,T}
  #
  approxConvOnElements!(ccwl, elements, ccwl.threadmodel)
end



"""
    $(SIGNATURES)

Prepare a common functor computation object `prepareCommonConvWrapper{T}` containing 
the user factor functor along with additional variables and information using during 
approximate convolution computations.

DevNotes
- TODO consolidate with others, see https://github.com/JuliaRobotics/IncrementalInference.jl/projects/6
"""
function prepareCommonConvWrapper!( F_::Type{<:AbstractRelative},
                                    ccwl::CommonConvWrapper{F},
                                    Xi::AbstractVector{<:DFGVariable},
                                    solvefor::Symbol,
                                    N::Int;
                                    needFreshMeasurements::Bool=true,
                                    solveKey::Symbol=:default  ) where {F <: AbstractFactor}
  #

  # FIXME, order of fmd ccwl cf are a little weird and should be revised.
  pttypes = getVariableType.(Xi) .|> getPointType
  PointType = 0 < length(pttypes) ? pttypes[1] : Vector{Float64}
  vecPtsArr = Vector{Vector{PointType}}()

  #TODO some better consolidate is needed
  ccwl.vartypes = typeof.(getVariableType.(Xi))

  # FIXME maxlen should parrot N (barring multi-/nullhypo issues)
  maxlen, sfidx, mani = prepareparamsarray!(vecPtsArr, Xi, solvefor, N, solveKey=solveKey)

  # FIXME ON FIRE, what happens if this is a partial dimension factor?  See #1246
  ccwl.xDim = getDimension(getVariableType(Xi[sfidx]))
  # ccwl.xDim = length(vecPtsArr[sfidx][1])
  # TODO should be selecting for the correct multihypothesis mode

    # setup the partial or complete decision variable dimensions for this ccwl object
    # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
    # TODO, should this not be part of `prepareCommonConvWrapper` -- only here do we look for .partial
    _setCCWDecisionDimsConv!(ccwl)

  # SHOULD WE SLICE ARR DOWN BY PARTIAL DIMS HERE (OR LATER)?
  ccwl.params = vecPtsArr # map( ar->view(ar, ccwl.partialDims, :), vecPtsArr)
  
  # get factor metadata -- TODO, populate, also see #784
  fmd = FactorMetadata(Xi, getLabel.(Xi), ccwl.params, solvefor, nothing)

  # TODO consolidate with ccwl??
  # FIXME do not divert Mixture for sampling
  # cf = _buildCalcFactorMixture(ccwl, fmd, 1, ccwl.measurement, ccwl.params) # TODO perhaps 0 is safer
  cf = CalcFactor( ccwl.usrfnc!, fmd, 0, length(ccwl.measurement), ccwl.measurement, ccwl.params)

  #  get variable node data
  vnds = Xi

  # option to disable fresh samples
  if needFreshMeasurements
    # TODO refactor
    ccwl.measurement = sampleFactor(cf, maxlen)
    # sampleFactor!(ccwl, maxlen, fmd, vnds)
  end

  if ccwl.specialzDim
    ccwl.zDim = ccwl.usrfnc!.zDim[sfidx]
  else
    ccwl.zDim = length(ccwl.measurement[1][1]) # FIXME -- zDim aspect needs to be reviewed
  end
  ccwl.varidx = sfidx

  # set each CPT
  for thrid in 1:Threads.nthreads()
    cpt_ = ccwl.cpt[thrid] 
    cpt_.X = ccwl.params[sfidx]

    # NOTE should be done during constructor for this factor only
    # resize!(cpt_.p, ccwl.xDim) 
    # cpt_.p .= 1:ccwl.xDim

    # used in ccw functor for AbstractRelativeMinimize
    # TODO JT - Confirm it should be updated here. Testing in prepgenericconvolution
    resize!(cpt_.res, ccwl.zDim) 
    fill!(cpt_.res, 0.0)
    # cpt_.res = zeros(ccwl.xDim) 
  end

  return sfidx, maxlen, mani
end


function prepareCommonConvWrapper!( ccwl::Union{CommonConvWrapper{F},
                                                CommonConvWrapper{Mixture{N_,F,S,T}}},
                                    Xi::AbstractVector{<:DFGVariable},
                                    solvefor::Symbol,
                                    N::Int;
                                    needFreshMeasurements::Bool=true,
                                    solveKey::Symbol=:default  ) where {N_,F<:AbstractRelative,S,T}
  #
  prepareCommonConvWrapper!(F, ccwl, Xi, solvefor, N, needFreshMeasurements=needFreshMeasurements, solveKey=solveKey)
end

function generateNullhypoEntropy( val::AbstractMatrix{<:Real},
                                  maxlen::Int,
                                  spreadfactor::Real=10  )
  #
  # covD = sqrt.(vec(Statistics.var(val,dims=2))) .+ 1e-3
  # cVar = diagm((spreadfactor*covD).^2)
  len = size(val, 1)
  cVar = diagm((spreadfactor*ones(len)).^2)
  mu = zeros( len )
  MvNormal( mu, cVar )
end

# FIXME SWITCH TO ON-MANIFOLD VERSION (see next, #TODO remove this one if manifold version is correct)
function calcVariableCovarianceBasic(M::AbstractManifold, ptsArr::Vector{Vector{Float64}})
  # cannot calculate the stdev from uninitialized state
  # FIXME assume point type is only Vector{Float} at this time
  @cast arr[i,j] := ptsArr[j][i]
  msst = Statistics.std(arr, dims=2)
  # FIXME use adaptive scale, see #802
  msst_ = 0 < sum(1e-10 .< msst) ? maximum(msst) : 1.0
  return msst_
end

function calcVariableCovarianceBasic(M::AbstractManifold, ptsArr::Vector{P}) where P
  #TODO double check the maths,. it looks like its working at least for groups
  μ = mean(M, ptsArr)
  Xcs = vee.(Ref(M), Ref(μ), log.(Ref(M), Ref(μ), ptsArr))
  Σ = mean(Xcs .* transpose.(Xcs))
  @debug "calcVariableCovarianceBasic" μ
  @debug "calcVariableCovarianceBasic" Σ
  # TODO don't know what to do here so keeping as before, #FIXME it will break
  # a change between this and previous is that a full covariance matrix is returned
  msst = Σ
  msst_ = 0 < sum(1e-10 .< msst) ? maximum(msst) : 1.0
  return msst_
end

"""
    $SIGNATURES

Control the amount of entropy to add to null-hypothesis in multihypo case.

Notes:
- Basically calculating the covariance (with a bunch of assumptions TODO, fix)
- FIXME, Currently only supports Euclidean domains.
- FIXME, allow particle subpopulations instead of just all of a variable
"""
function calcVariableDistanceExpectedFractional(ccwl::CommonConvWrapper,
                                                sfidx::Int,
                                                certainidx::Vector{Int};
                                                kappa::Float64=3.0  )
  #
  if sfidx in certainidx
    msst_ = calcVariableCovarianceBasic(getManifold(ccwl.vartypes[sfidx]), ccwl.params[sfidx])
    return kappa*msst_
  end
  # @assert !(sfidx in certainidx) "null hypo distance does not work for sfidx in certainidx"

  # get mean of all fractional variables
  # ccwl.params::Vector{Vector{P}}
  uncertainidx = setdiff(1:length(ccwl.params), certainidx)
  uncMeans = zeros(length(ccwl.params[sfidx][1]), length(uncertainidx))
  dists = zeros(length(uncertainidx)+length(certainidx))
  dims = length(ccwl.params[sfidx][1])
  for (count,i) in enumerate(uncertainidx)
    # uncMeans[:,count] = Statistics.mean(ccwl.params[i], dims=2)[:]
    for pr in ccwl.params[i]
      uncMeans[:,count] .+= pr
    end
    uncMeans[:,count] ./= length(ccwl.params[i])
  end
  count = 0
  # refMean = Statistics.mean(ccwl.params[sfidx], dims=2)[:]
  refMean = zeros(length(ccwl.params[sfidx][1]))
  for pr in ccwl.params[sfidx]
    refMean .+= pr
  end
  refMean ./= length(ccwl.params[sfidx])

  # calc for uncertain and certain
  for i in uncertainidx
    count += 1
    dists[count] = norm(refMean - uncMeans[:,count])
  end
  # also check distance to certainidx for general scale reference (workaround heuristic)
  for cidx in certainidx
    count += 1
    cerMean = zeros(length(ccwl.params[cidx][1]))
    # cerMean = Statistics.mean(ccwl.params[cidx], dims=2)[:]
    for pr in ccwl.params[cidx]
      cerMean .+= pr
    end
    cerMean ./= length(ccwl.params[cidx])
    dists[count] = norm(refMean[1:dims] - cerMean[1:dims])
  end

  push!(dists, 1e-2)
  return kappa*maximum(dists)
end

# Add entrypy on a point in `points` on manifold M, only on dimIdx if in p 
function addEntropyOnManifoldHack!( M::ManifoldsBase.AbstractManifold,
                                    points::Union{<:AbstractVector{<:Real},SubArray}, 
                                    dimIdx::AbstractVector, 
                                    spreadDist::Real,
                                    p::Union{Colon, <:AbstractVector}=: )
  #
  if length(points) == 0 
    return nothing
  end
  
  # preallocate 
  T = number_eltype(points[1])
  Xc = zeros(T, manifold_dimension(M))
  X = get_vector(M, points[1], Xc, DefaultOrthogonalBasis())  
  
  for idx in 1:length(points)
    # build tangent coordinate random
    for dim in dimIdx
      if (p === :) || dim in p
        Xc[dim] = spreadDist*(rand()-0.5)
      end
    end
    # update tangent vector X
    get_vector!(M, X, points[idx], Xc, DefaultOrthogonalBasis())  
    #update point
    exp!(M, points[idx], points[idx], X)
    
  end
  #
  # manis = convert(Tuple, M) # LEGACY, TODO REMOVE
  # # TODO deprecate
  # maniAddOps, _, _, _ = buildHybridManifoldCallbacks(manis)
  # # add 1σ "noise" level to max distance as control
  # # 1:size(addEntr, 1)
  # for dim in dimIdx, idx in 1:length(addEntr)
  #   if (p === :) || dim in p
  #     addEntr[idx][dim] = maniAddOps[dim](addEntr[idx][dim], spreadDist*(rand()-0.5))
  #   end
  # end
  nothing
end

"""
    $(SIGNATURES)

Common function to compute across a single user defined multi-hypothesis ambiguity per factor.  
This function dispatches both `AbstractRelativeRoots` and `AbstractRelativeMinimize` factors.
"""
function computeAcrossHypothesis!(ccwl::Union{<:CommonConvWrapper{F},
                                              <:CommonConvWrapper{Mixture{N_,F,S,T}}},
                                  allelements::AbstractVector,
                                  activehypo,
                                  certainidx::Vector{Int},
                                  sfidx::Int,
                                  maxlen::Int,
                                  mani::ManifoldsBase.AbstractManifold; # maniAddOps::Tuple;
                                  spreadNH::Real=5.0,
                                  inflateCycles::Int=3,
                                  skipSolve::Bool=false ) where {N_,F<:AbstractRelative,S,T}
  #
  count = 0

  cpt_ = ccwl.cpt[Threads.threadid()]
  
  # @assert norm(ccwl.certainhypo - certainidx) < 1e-6
  for (hypoidx, vars) in activehypo
    count += 1
    
    # now do hypothesis specific
    if sfidx in certainidx && hypoidx != 0 || hypoidx in certainidx || hypoidx == sfidx
      # hypo case hypoidx, sfidx = $hypoidx, $sfidx
      for i in 1:Threads.nthreads()  ccwl.cpt[i].activehypo = vars; end
      
      addEntr = view(ccwl.params[sfidx], allelements[count])
      # dynamic estimate with user requested speadNH of how much noise to inject (inflation or nullhypo)
      spreadDist = calcVariableDistanceExpectedFractional(ccwl, sfidx, certainidx, kappa=ccwl.inflation)
      
      # do proposal inflation step, see #1051
      # consider duplicate convolution approximations for inflation off-zero
      # ultimately set by dfg.params.inflateCycles
      for iflc in 1:inflateCycles
        addEntropyOnManifoldHack!(mani, addEntr, 1:getDimension(mani), spreadDist, cpt_.p)
        # no calculate new proposal belief on kernels `allelements[count]`
        skipSolve ? @warn("skipping numerical solve operation") : approxConvOnElements!(ccwl, allelements[count])
      end
    elseif hypoidx != sfidx && hypoidx != 0
      # snap together case
      # multihypo, take other value case
      # sfidx=2, hypoidx=3:  2 should take a value from 3
      # sfidx=3, hypoidx=2:  3 should take a value from 2
      # DEBUG sfidx=2, hypoidx=1 -- bad when do something like multihypo=[0.5;0.5] -- issue 424
      # ccwl.params[sfidx][:,allelements[count]] = view(ccwl.params[hypoidx],:,allelements[count])
        # NOTE make alternative case only operate as null hypo
        addEntr = view(ccwl.params[sfidx], allelements[count])
        # dynamic estimate with user requested speadNH of how much noise to inject (inflation or nullhypo)
        spreadDist = calcVariableDistanceExpectedFractional(ccwl, sfidx, certainidx, kappa=spreadNH)
        addEntropyOnManifoldHack!(mani, addEntr, 1:getDimension(mani), spreadDist)

    elseif hypoidx == 0
      # basically do nothing since the factor is not active for these allelements[count]
      # inject more entropy in nullhypo case
      # add noise (entropy) to spread out search in convolution proposals
      addEntr = view(ccwl.params[sfidx], allelements[count])
      # dynamic estimate with user requested speadNH of how much noise to inject (inflation or nullhypo)
      spreadDist = calcVariableDistanceExpectedFractional(ccwl, sfidx, certainidx, kappa=spreadNH)
      # # make spread (1σ) equal to mean distance of other fractionals
      addEntropyOnManifoldHack!(mani, addEntr, 1:getDimension(mani), spreadDist)
    else
      error("computeAcrossHypothesis -- not dealing with multi-hypothesis case correctly")
    end
  end
  nothing
end
  # elseif hypoidx == sfidx
  #   # multihypo, do conv case, hypoidx == sfidx
  #   ah = sort(union([sfidx;], certainidx))
  #   @assert norm(ah - vars) < 1e-10
  #   for i in 1:Threads.nthreads()  ccwl.cpt[i].activehypo = ah; end
  #   approxConvOnElements!(ccwl, allelements[count])


"""
    $SIGNATURES
Internal method to set which dimensions should be used as the decision variables for later numerical optimization.
"""
function _setCCWDecisionDimsConv!(ccwl::Union{CommonConvWrapper{F},
                                              CommonConvWrapper{Mixture{N_,F,S,T}}} ) where {N_,F<:Union{AbstractManifoldMinimize, AbstractRelativeMinimize, AbstractRelativeRoots, AbstractPrior},S,T}
  #
  # return nothing

  p = if ccwl.partial
    Int32[ccwl.usrfnc!.partial...]
  else
    Int32[1:ccwl.xDim...]
  end

  ccwl.partialDims = (p)
  # NOTE should only be done in the constructor
  for thrid in 1:Threads.nthreads()
    length(ccwl.cpt[thrid].p) != length(p) ? resize!(ccwl.cpt[thrid].p, length(p)) : nothing
    ccwl.cpt[thrid].p .= p # SVector... , see ccw.partialDims
  end
  nothing
end

"""
    $(SIGNATURES)

Multiple dispatch wrapper for `<:AbstractRelativeRoots` types, to prepare and execute the general approximate convolution with user defined factor residual functions.  This method also supports multihypothesis operations as one mechanism to introduce new modality into the proposal beliefs.

Planned changes will fold null hypothesis in as a standard feature and no longer appear as a separate `InferenceType`.
"""
function evalPotentialSpecific( Xi::AbstractVector{<:DFGVariable},
                                ccwl::CommonConvWrapper{T},
                                solvefor::Symbol,
                                T_::Type{<:AbstractRelative},
                                measurement::Tuple=(zeros(0,100),);
                                needFreshMeasurements::Bool=true,
                                solveKey::Symbol=:default,
                                N::Int=size(measurement[1],2),
                                spreadNH::Real=3.0,
                                inflateCycles::Int=3,
                                dbg::Bool=false,
                                skipSolve::Bool=false  ) where {T <: AbstractFactor}
  #

  # Prep computation variables
  # NOTE #1025, should FMD be built here...
  sfidx, maxlen, mani = prepareCommonConvWrapper!(ccwl, Xi, solvefor, N, needFreshMeasurements=needFreshMeasurements, solveKey=solveKey)
  # check for user desired measurement values
  if 0 < length(measurement[1])
    ccwl.measurement = measurement
  end

  # Check which variables have been initialized
  isinit = map(x->isInitialized(x), Xi)
  
  
  # assemble how hypotheses should be computed
  # TODO convert to HypothesisRecipeElements result
  _, allelements, activehypo, mhidx = assembleHypothesesElements!(ccwl.hypotheses, maxlen, sfidx, length(Xi), isinit, ccwl.nullhypo )
  certainidx = ccwl.certainhypo
  
  # get manifold add operations
  # TODO, make better use of dispatch, see JuliaRobotics/RoME.jl#244
  # addOps, d1, d2, d3 = buildHybridManifoldCallbacks(manis)
  mani = getManifold(getVariableType(Xi[sfidx]))

  # perform the numeric solutions on the indicated elements
  # error("ccwl.xDim=$(ccwl.xDim)")
  # FIXME consider repeat solve as workaround for inflation off-zero 
  computeAcrossHypothesis!( ccwl, allelements, activehypo, certainidx, 
                            sfidx, maxlen, mani, spreadNH=spreadNH, 
                            inflateCycles=inflateCycles, skipSolve=skipSolve )
  #
  return ccwl.params[ccwl.varidx]
end


# TODO `measurement` might not be properly wired up yet
# TODO consider 1051 here to inflate proposals as general behaviour
function evalPotentialSpecific( Xi::AbstractVector{<:DFGVariable},
                                ccwl::CommonConvWrapper{T},
                                solvefor::Symbol,
                                T_::Type{<:AbstractPrior},
                                measurement::Tuple=(Vector{Vector{Float64}}(),);
                                needFreshMeasurements::Bool=true,
                                solveKey::Symbol=:default,
                                N::Int=length(measurement[1]),
                                dbg::Bool=false,
                                spreadNH::Real=3.0,
                                inflateCycles::Int=3,
                                skipSolve::Bool=false ) where {T <: AbstractFactor}
  #
  # setup the partial or complete decision variable dimensions for this ccwl object
  # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
  _setCCWDecisionDimsConv!(ccwl)
  
  # FIXME, NEEDS TO BE CLEANED UP AND WORK ON MANIFOLDS PROPER
  fnc = ccwl.usrfnc!
  sfidx = findfirst(getLabel.(Xi) .== solvefor)
  # sfidx = 1 #  WHY HARDCODED TO 1??
  solveForPts = getVal(Xi[sfidx], solveKey=solveKey)
  nn = maximum([N; length(measurement[1]); length(solveForPts); length(ccwl.params[sfidx])])
  # vnds = Xi # (x->getSolverData(x)).(Xi)
  # FIXME better standardize in-place operations (considering solveKey)
  if needFreshMeasurements
    cf = CalcFactor( ccwl )
    ccwl.measurement = sampleFactor(cf, nn)
  end
  # Check which variables have been initialized
  isinit::Vector{Bool} = Xi .|> isInitialized .|> Bool
  _, _, _, mhidx = assembleHypothesesElements!(ccwl.hypotheses, nn, sfidx, length(Xi), isinit, ccwl.nullhypo )
  # get solvefor manifolds, FIXME ON FIRE, upgrade to new Manifolds.jl
  mani = getManifold(Xi[sfidx])
  # two cases on how to use the measurement
  nhmask = mhidx .== 0
  ahmask = mhidx .== 1
  # generate nullhypo samples
  # inject lots of entropy in nullhypo case
  # make spread (1σ) equal to mean distance of other fractionals
  # FIXME better standardize in-place operations (considering solveKey)
  addEntr = if length(solveForPts) == nn
    deepcopy(solveForPts)
  else
    ret = typeof(solveForPts)(undef, nn)
    for i in 1:length(solveForPts)
      ret[i] = solveForPts[i]
    end
    for i in (length(solveForPts)+1):nn
      ret[i] = getPointIdentity(getVariableType(Xi[sfidx]))
    end
    ret
  end
  # @show nn, size(addEntr), size(nhmask), size(solveForPts)
  addEntrNH = view(addEntr, nhmask)
  spreadDist = spreadNH*calcVariableCovarianceBasic(mani, addEntr)
  # partials are treated differently
  if !ccwl.partial
      # TODO for now require measurements to be coordinates too
      # @show typeof(ccwl.measurement[1])
      for m in (1:length(ahmask))[ahmask]
        # @show m, addEntr[m], ccwl.measurement[1][m]
        addEntr[m] .= ccwl.measurement[1][m]
      end
      # ongoing part of RoME.jl #244
      addEntropyOnManifoldHack!(mani, addEntrNH, 1:getDimension(mani), spreadDist)
  else
    i = 0
    for dimnum in fnc.partial
      i += 1
      for m in (1:length(ahmask))[ahmask]
        addEntr[m][dimnum] = ccwl.measurement[1][m][i]
      end
      # @show size(addEntr), dimnum, nhmask
      addEntrNHp = view(view(addEntr, (1:length(ahmask))[ahmask]), dimnum)
      # ongoing part of RoME.jl #244
      addEntropyOnManifoldHack!(mani, addEntrNHp, dimnum:dimnum, spreadDist)
    end
  end
  return addEntr
end


function evalPotentialSpecific( Xi::AbstractVector{<:DFGVariable},
                                ccwl::CommonConvWrapper{Mixture{N_,F,S,T}},
                                solvefor::Symbol,
                                measurement::Tuple=(Vector{Vector{Float64}}(),);
                                needFreshMeasurements::Bool=true,
                                solveKey::Symbol=:default,
                                N::Int=length(measurement[1]),
                                dbg::Bool=false,
                                spreadNH::Real=3.0,
                                inflateCycles::Int=3,
                                skipSolve::Bool=false ) where {N_,F<:AbstractFactor,S,T}
  #
  evalPotentialSpecific(Xi,
                        ccwl,
                        solvefor,
                        F,
                        measurement;
                        needFreshMeasurements=needFreshMeasurements,
                        solveKey=solveKey,
                        N=N,
                        dbg=dbg,
                        spreadNH=spreadNH,
                        inflateCycles=inflateCycles,
                        skipSolve=skipSolve )
end


function evalPotentialSpecific( Xi::AbstractVector{<:DFGVariable},
                                ccwl::CommonConvWrapper{F},
                                solvefor::Symbol,
                                measurement::Tuple=(Vector{Vector{Float64}}(),);
                                needFreshMeasurements::Bool=true,
                                solveKey::Symbol=:default,
                                N::Int=length(measurement[1]),
                                dbg::Bool=false,
                                spreadNH::Real=3.0,
                                inflateCycles::Int=3,
                                skipSolve::Bool=false ) where {F <: AbstractFactor}
  #
  evalPotentialSpecific(Xi,
                        ccwl,
                        solvefor,
                        F,
                        measurement;
                        needFreshMeasurements=needFreshMeasurements,
                        solveKey=solveKey,
                        N=N,
                        dbg=dbg,
                        spreadNH=spreadNH,
                        inflateCycles=inflateCycles,
                        skipSolve=skipSolve )
end


"""
    $(SIGNATURES)

Single entry point for evaluating factors from factor graph, using multiple dispatch to locate the correct `evalPotentialSpecific` function.
"""
function evalFactor(dfg::AbstractDFG,
                    fct::DFGFactor,
                    solvefor::Symbol,
                    measurement::Tuple=(Vector{Vector{Float64}}(),);
                    needFreshMeasurements::Bool=true,
                    solveKey::Symbol=:default,
                    N::Int=length(measurement[1]),
                    inflateCycles::Int=getSolverParams(dfg).inflateCycles,
                    dbg::Bool=false,
                    skipSolve::Bool=false  )
  #

  ccw = _getCCW(fct)
  # TODO -- this build up of Xi is excessive and could happen at addFactor time
  variablelist = getVariableOrder(fct)
  Xi = getVariable.(dfg, variablelist)

  # setup operational values before compute (likely to be refactored) 
  for i in 1:Threads.nthreads()
    ccw.cpt[i].factormetadata.variablelist = variablelist
    ccw.cpt[i].factormetadata.solvefor = solvefor
  end

  return evalPotentialSpecific( Xi, ccw, solvefor, measurement, needFreshMeasurements=needFreshMeasurements,
                                solveKey=solveKey, N=N, dbg=dbg, spreadNH=getSolverParams(dfg).spreadNH, 
                                inflateCycles=inflateCycles, skipSolve=skipSolve )
  #
end

"""
    $SIGNATURES

Helper function for evaluating factor residual functions, by adding necessary `CalcFactor` wrapper.
  
Notes
- Factor must already be in a factor graph to work
- Will not yet properly support all multihypo nuances, more a function for testing
- Useful for debugging a factor. 

Example
```julia
fg = generateCanonicalFG_Kaess()

residual = calcFactorResidual(fg, :x1x2f1, [1.0], [0.0], [0.0])
```
"""
calcFactorResidual(dfg::AbstractDFG, fctsym::Symbol, args...) = CalcFactor(IIF._getCCW(dfg, fctsym))(args...)


function approxConv(dfg::AbstractDFG,
                    fc::DFGFactor,
                    target::Symbol,
                    measurement::Tuple=(Vector{Vector{Float64}}(),);
                    solveKey::Symbol=:default,
                    N::Int=length(measurement[1]), 
                    skipSolve::Bool=false )
  #
  v1 = getVariable(dfg, target)
  N = N == 0 ? getNumPts(v1) : N
  return evalFactor(dfg, fc, v1.label, measurement, solveKey=solveKey, N=N, skipSolve=skipSolve)
end


"""
    $SIGNATURES

Calculate the sequential series of convolutions in order as listed by `fctLabels`, and starting from the 
value already contained in the first variable.  

Notes
- `target` must be a variable.
- The ultimate `target` variable must be given to allow path discovery through n-ary factors.
- Fresh starting point will be used if first element in `fctLabels` is a unary `<:AbstractPrior`.
- This function will not change any values in `dfg`, and might have slightly less speed performance to meet this requirement.
- pass in `tfg` to get a recoverable result of all convolutions in the chain.
- `setPPE` and `setPPEmethod` can be used to store PPE information in temporary `tfg`

DevNotes
- TODO strong requirement that this function is super efficient on single factor/variable case!
- FIXME must consolidate with `accumulateFactorMeans`
- TODO `solveKey` not fully wired up everywhere yet
  - tfg gets all the solveKeys inside the source `dfg` variables
- TODO add a approxConv on PPE option
  - Consolidate with [`accumulateFactorMeans`](@ref), `approxConvBinary`

Related

[`approxDeconv`](@ref), `LightDFG.findShortestPathDijkstra`, [`evalFactor`](@ref)
"""
function approxConv(dfg::AbstractDFG, 
                    from::Symbol, 
                    target::Symbol,
                    measurement::Tuple=(Vector{Vector{Float64}}(),);
                    solveKey::Symbol=:default,
                    N::Int = length(measurement[1]),
                    tfg::AbstractDFG = initfg(),
                    setPPEmethod::Union{Nothing, Type{<:AbstractPointParametricEst}}=nothing,
                    setPPE::Bool= setPPEmethod !== nothing,
                    path::AbstractVector{Symbol}=Symbol[],
                    skipSolve::Bool=false  )
  #
  # @assert isVariable(dfg, target) "approxConv(dfg, from, target,...) where `target`=$target must be a variable in `dfg`"
  
  if from in ls(dfg, target)
    # direct request
    # TODO avoid this allocation for direct cases ( dfg, :x1x2f1, :x2[/:x1] )
    path = Symbol[from; target]
    varLbls = Symbol[target;]
  else
    # must first discover shortest factor path in dfg
    # TODO DFG only supports LightDFG.findShortestPathDijkstra at the time of writing (DFG v0.10.9)
    path = 0 == length(path) ? findShortestPathDijkstra(dfg, from, target) : path
    @assert path[1] == from "sanity check that shortest path function is working as expected"

    # list of variables
    fctMsk = isFactor.(dfg, path)
    # which factors in the path
    fctLbls = path[fctMsk]
    # must still add
    varLbls =  union(lsf.(dfg, fctLbls)...)
    neMsk = exists.(tfg, varLbls) .|> x-> xor(x,true)
    # put the non-existing variables into the temporary graph `tfg`
    # bring all the solveKeys too
    addVariable!.(tfg, getVariable.(dfg, varLbls[neMsk]))
    # variables adjacent to the shortest path should be initialized from dfg
    setdiff(varLbls, path[xor.(fctMsk,true)]) .|> x->initManual!(tfg, x, getBelief(dfg, x))
  end
  
  # find/set the starting point
  idxS = 1
  pts = if varLbls[1] == from
    # starting from a variable
    pts0 = getBelief(dfg, varLbls[1]) |> getPoints
  else
    # chain would start one later
    idxS += 1
    # get the factor
    fct0 = getFactor(dfg,from)
    # get the Matrix{<:Real} of projected points
    pts1 = approxConv(dfg, fct0, path[2], measurement, solveKey=solveKey, N=N, skipSolve=skipSolve)
    length(path) == 2 ? (return pts1) : pts1
  end
  # didn't return early so shift focus to using `tfg` more intensely
  initManual!(tfg, varLbls[1], pts)
  # use in combination with setPPE and setPPEmethod keyword arguments
  ppemethod = setPPEmethod === nothing ? MeanMaxPPE : setPPEmethod
  !setPPE ? nothing : setPPE!(tfg, varLbls[1], solveKey, ppemethod)

  # do chain of convolutions
  for idx in idxS:length(path)
    if fctMsk[idx]
      # this is a factor path[idx]
      fct = getFactor(dfg, path[idx])
      addFactor!(tfg, fct)
      pts = approxConv(tfg, fct, path[idx+1], solveKey=solveKey, N=N, skipSolve=skipSolve)
      initManual!(tfg, path[idx+1], pts)
      !setPPE ? nothing : setPPE!(tfg, path[idx+1], solveKey, ppemethod)
    end
  end

  # return target variable values
  return getBelief(tfg, target) |> getPoints
end



## ====================================================================================
## TODO better consolidate below with existing functions
## ====================================================================================




# TODO should this be consolidated with regular approxConv?
# TODO, perhaps pass Xi::Vector{DFGVariable} instead?
function approxConvBinary(arr::Vector{Vector{Float64}},
                          meas::AbstractFactor,
                          outdims::Int,
                          fmd::FactorMetadata,
                          measurement::Tuple=(Vector{Vector{Float64}}(),);
                          varidx::Int=2,
                          N::Int=length(arr),
                          vnds=DFGVariable[] )
  #
  # N = N == 0 ? size(arr,2) : N
  pts = [zeros(outdims) for _ in 1:N];
  ptsArr = Vector{Vector{Vector{Float64}}}()
  push!(ptsArr,arr)
  push!(ptsArr,pts)

  fmd.arrRef = ptsArr

  # TODO consolidate with ccwl??
  # FIXME do not divert Mixture for sampling
  # cf = _buildCalcFactorMixture(ccwl, fmd, 1, ccwl.measurement, ARR) # TODO perhaps 0 is safer
  # FIXME 0, 0, ()
  cf = CalcFactor( meas, fmd, 0, 0, (), ptsArr)

  measurement = length(measurement[1]) == 0 ? sampleFactor(cf, N) : measurement
  # measurement = size(measurement[1],2) == 0 ? sampleFactor(meas, N, fmd, vnds) : measurement

  zDim = length(measurement[1][1])
  ccw = CommonConvWrapper(meas, ptsArr[varidx], zDim, ptsArr, fmd, varidx=varidx, measurement=measurement)  # N=> size(measurement[1],2)

  for n in 1:N
    ccw.cpt[Threads.threadid()].particleidx = n
    _solveCCWNumeric!( ccw )
  end
  return pts
end



"""
    $SIGNATURES

Calculate both measured and predicted relative variable values, starting with `from` at zeros up to `to::Symbol`.

Notes
- assume single variable separators only.
"""
function accumulateFactorChain( dfg::AbstractDFG,
                                from::Symbol,
                                to::Symbol,
                                fsyms::Vector{Symbol}=findFactorsBetweenNaive(dfg, from, to);
                                initval=zeros(size(getVal(dfg, from))))

  # get associated variables
  svars = union(ls.(dfg, fsyms)...)

  # use subgraph copys to do calculations
  tfg_meas = buildSubgraph(dfg, [svars;fsyms])
  tfg_pred = buildSubgraph(dfg, [svars;fsyms])

  # drive variable values manually to ensure no additional stochastics are introduced.
  nextvar = from
  initManual!(tfg_meas, nextvar, initval)
  initManual!(tfg_pred, nextvar, initval)

  # nextfct = fsyms[1] # for debugging
  for nextfct in fsyms
    nextvars = setdiff(ls(tfg_meas,nextfct),[nextvar])
    @assert length(nextvars) == 1 "accumulateFactorChain requires each factor pair to separated by a single variable"
    nextvar = nextvars[1]
    meas, pred = approxDeconv(dfg, nextfct) # solveFactorMeasurements
    pts_meas = approxConv(tfg_meas, nextfct, nextvar, (meas,ones(Int,100),collect(1:100)))
    pts_pred = approxConv(tfg_pred, nextfct, nextvar, (pred,ones(Int,100),collect(1:100)))
    initManual!(tfg_meas, nextvar, pts_meas)
    initManual!(tfg_pred, nextvar, pts_pred)
  end
  return getVal(tfg_meas,nextvar), getVal(tfg_pred,nextvar)
end




"""
    $(SIGNATURES)

Compute proposal belief on `vertid` through `fct` representing some constraint in factor graph.
Always full dimension variable node -- partial constraints will only influence subset of variable dimensions.
The remaining dimensions will keep pre-existing variable values.

Notes
- fulldim is true when "rank-deficient" -- TODO swap to false (or even float)
"""
function findRelatedFromPotential(dfg::AbstractDFG,
                                  fct::DFGFactor,
                                  target::Symbol,
                                  measurement::Tuple=(Vector{Vector{Float64}}(),);
                                  N::Int=length(measurement[1]),
                                  solveKey::Symbol=:default,
                                  dbg::Bool=false  )
  #

  # # assuming it is properly initialized TODO
  pts = evalFactor(dfg, fct, target, solveKey=solveKey, N=N, dbg=dbg);
  # pts = approxConv(dfg, fct, target, measurement, N=N, solveKey=solveKey)
  
  # # determine if evaluation is "dimension-deficient"
  # solvable dimension
  inferdim = getFactorSolvableDim(dfg, fct, target, solveKey)
  # zdim = getFactorDim(fct)
  # vdim = getVariableDim(DFG.getVariable(dfg, target))

  # TODO -- better to upsample before the projection
  # Ndim = length(pts[1]) # not used, should maybe be: 
  # Ndim = manifold_dimension(M)
  Npoints = length(pts) # size(pts,2)
  # Assume we only have large particle population sizes, thanks to addNode!
  M = getManifold(getVariableType(dfg, target))
  proposal = AMP.manikde!(M, pts)

  # FIXME consolidate with approxConv method instead
  if Npoints != N # this is where we control the overall particle set size
      proposal = AMP.resample(proposal,N)
  end
  return (proposal, inferdim)
end


# function _expandNamedTupleType()
  
# end


"""
    $SIGNATURES

Compute the proposals of a destination vertex for each of `factors` and place the result
as belief estimates in both `dens` and `partials` respectively.

Notes
- TODO: also return if proposals were "dimension-deficient" (aka ~rank-deficient).
"""
function proposalbeliefs!(dfg::AbstractDFG,
                          destvertlabel::Symbol,
                          factors::AbstractVector{<:DFGFactor},
                          dens::Vector{<:ManifoldKernelDensity},
                          # partials::Dict{Any, Vector{ManifoldKernelDensity}}, # TODO change this structure
                          measurement::Tuple=(Vector{Vector{Float64}}(),);
                          solveKey::Symbol=:default,
                          N::Int=100,
                          dbg::Bool=false  )
  #
  

  # group partial dimension factors by selected dimensions -- i.e. [(1,)], [(1,2),(1,2)], [(2,);(2;)]


  # populate the full and partial dim containers
  inferddimproposal = Vector{Float64}(undef, length(factors))
  for (count,fct) in enumerate(factors)
    data = getSolverData(fct)
    ccwl = _getCCW(data)
    propBel_, inferd = findRelatedFromPotential(dfg, fct, destvertlabel, measurement, N=N, dbg=dbg, solveKey=solveKey)
    # partial density
    propBel = if isPartial(ccwl)
      pardims = _getDimensionsPartial(ccwl)
      @assert [getFactorType(fct).partial...] == [pardims...] "partial dims error $(getFactorType(fct).partial) vs $pardims"
      AMP.marginal(propBel_, Int[pardims...])
    else
      propBel_
    end
    push!(dens, propBel)
    inferddimproposal[count] = inferd
  end
  inferddimproposal
end



#
