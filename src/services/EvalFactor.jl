

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
                                elements::Union{Vector{Int}, UnitRange{Int}}, 
                                ::Type{<:MultiThreaded},
                                _slack=nothing ) where {N_,F<:AbstractRelative,S,T}
  #
  Threads.@threads for n in elements
    # ccwl.thrid_ = Threads.threadid()
    ccwl.cpt[Threads.threadid()].particleidx = n
    
    # ccall(:jl_, Nothing, (Any,), "starting loop, thrid_=$(Threads.threadid()), partidx=$(ccwl.cpt[Threads.threadid()].particleidx)")
    _solveCCWNumeric!( ccwl, _slack=_slack)
  end
  nothing
end


function approxConvOnElements!( ccwl::Union{CommonConvWrapper{F},
                                            CommonConvWrapper{Mixture{N_,F,S,T}}},
                                elements::Union{Vector{Int}, UnitRange{Int}}, 
                                ::Type{<:SingleThreaded},
                                _slack=nothing ) where {N_,F<:AbstractRelative,S,T}
  #
  for n in elements
    ccwl.cpt[Threads.threadid()].particleidx = n
    _solveCCWNumeric!( ccwl, _slack=_slack)
  end
  nothing
end


function approxConvOnElements!( ccwl::Union{CommonConvWrapper{F},
                                            CommonConvWrapper{Mixture{N_,F,S,T}}},
                                elements::Union{Vector{Int}, UnitRange{Int}},
                                _slack=nothing )  where {N_,F<:AbstractRelative,S,T}
  #
  approxConvOnElements!(ccwl, elements, ccwl.threadmodel, _slack)
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
    msst_ = sqrt(calcCovarianceBasic(getManifold(ccwl.vartypes[sfidx]), ccwl.params[sfidx]))
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
function addEntropyOnManifold!( M::ManifoldsBase.AbstractManifold,
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
                                  skipSolve::Bool=false,
                                  _slack=nothing ) where {N_,F<:AbstractRelative,S,T}
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
        addEntropyOnManifold!(mani, addEntr, 1:getDimension(mani), spreadDist, cpt_.p)
        # no calculate new proposal belief on kernels `allelements[count]`
        skipSolve ? @warn("skipping numerical solve operation") : approxConvOnElements!(ccwl, allelements[count], _slack)
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
        addEntropyOnManifold!(mani, addEntr, 1:getDimension(mani), spreadDist)

    elseif hypoidx == 0
      # basically do nothing since the factor is not active for these allelements[count]
      # inject more entropy in nullhypo case
      # add noise (entropy) to spread out search in convolution proposals
      addEntr = view(ccwl.params[sfidx], allelements[count])
      # dynamic estimate with user requested speadNH of how much noise to inject (inflation or nullhypo)
      spreadDist = calcVariableDistanceExpectedFractional(ccwl, sfidx, certainidx, kappa=spreadNH)
      # # make spread (1σ) equal to mean distance of other fractionals
      addEntropyOnManifold!(mani, addEntr, 1:getDimension(mani), spreadDist)
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
    $(SIGNATURES)

Multiple dispatch wrapper for `<:AbstractRelativeRoots` types, to prepare and execute the general approximate convolution with user defined factor residual functions.  This method also supports multihypothesis operations as one mechanism to introduce new modality into the proposal beliefs.

Planned changes will fold null hypothesis in as a standard feature and no longer appear as a separate `InferenceType`.
"""
function evalPotentialSpecific( Xi::AbstractVector{<:DFGVariable},
                                ccwl::CommonConvWrapper{T},
                                solvefor::Symbol,
                                T_::Type{<:AbstractRelative},
                                measurement::Tuple=(Vector{Float64}(),);
                                needFreshMeasurements::Bool=true,
                                solveKey::Symbol=:default,
                                N::Int= 0<length(measurement[1]) ? length(measurement[1]) : maximum(Npts.(getBelief.(Xi, solveKey))),
                                spreadNH::Real=3.0,
                                inflateCycles::Int=3,
                                dbg::Bool=false,
                                skipSolve::Bool=false,
                                _slack=nothing  ) where {T <: AbstractFactor}
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
  # FIXME consider repeat solve as workaround for inflation off-zero 
  computeAcrossHypothesis!( ccwl, allelements, activehypo, certainidx, 
                            sfidx, maxlen, mani, spreadNH=spreadNH, 
                            inflateCycles=inflateCycles, skipSolve=skipSolve,
                            _slack=_slack )
  #
  # do info per coord
  ipc = if ccwl._gradients === nothing 
    ones(getDimension(Xi[sfidx]))
  else
    ipc_ = ones(getDimension(Xi[sfidx]))
    # calcPerturbationFromVariable(ccwl._gradients, 2, ipc_) # TODO, WIP
    ipc_                                                     # TODO, WIP
  end

  # return the found points, and info per coord
  return ccwl.params[ccwl.varidx], ipc
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
                                skipSolve::Bool=false,
                                _slack=nothing ) where {T <: AbstractFactor}
  #
  # setup the partial or complete decision variable dimensions for this ccwl object
  # NOTE perhaps deconv has changed the decision variable list, so placed here during consolidation phase
  _setCCWDecisionDimsConv!(ccwl)

  # FIXME, NEEDS TO BE CLEANED UP AND WORK ON MANIFOLDS PROPER
  fnc = ccwl.usrfnc!
  sfidx = findfirst(getLabel.(Xi) .== solvefor)
  # sfidx = 1 #  WHY HARDCODED TO 1??
  solveForPts = getVal(Xi[sfidx], solveKey=solveKey)
  nn = maximum([N; calcZDim(ccwl); length(solveForPts); length(ccwl.params[sfidx])])  # length(measurement[1])

  # FIXME better standardize in-place operations (considering solveKey)
  if needFreshMeasurements
    cf = CalcFactor( ccwl )
    newMeas = sampleFactor(cf, nn)
    ccwl.measurement = newMeas
  end

  # Check which variables have been initialized, TODO not sure why forcing to Bool vs BitVector
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

  # view on elements marked for nullhypo
  addEntrNH = view(addEntr, nhmask)
  spreadDist = spreadNH*sqrt(calcCovarianceBasic(mani, addEntr))
  # partials are treated differently
  ipc = if !isPartial(ccwl) #ccwl.partial
    # TODO for now require measurements to be coordinates too
    # @show typeof(ccwl.measurement[1])
    for m in (1:length(addEntr))[ahmask]
      # FIXME, selection for all measurement::Tuple elements
      # @info "check broadcast" ccwl.usrfnc! addEntr[m] ccwl.measurement[1][m]
      setPointsMani!(addEntr[m], ccwl.measurement[1][m])
    end
    # ongoing part of RoME.jl #244
    addEntropyOnManifold!(mani, addEntrNH, 1:getDimension(mani), spreadDist)
    # do info per coords
    ones(getDimension(Xi[sfidx]))
  else
    # FIXME but how to add partial factor info only on affected dimensions fro general manifold points?
    pvec = [fnc.partial...]
    # active hypo that receives the regular measurement information
    for m in (1:length(addEntr))[ahmask]
      # addEntr is no longer in coordinates, these are now general manifold points!!
      # for (i,dimnum) in enumerate(fnc.partial)
        # FIXME, need ability to replace partial points
        partialCoords = ccwl.partialDims
        Msrc, = getManifoldPartial(mani,partialCoords)
        setPointsManiPartial!(mani, addEntr[m], Msrc, ccwl.measurement[1][m], partialCoords)
        # addEntr[m][dimnum] = ccwl.measurement[1][m][i]
      # end
    end
    # null hypo mask that needs to be perturbed by "noise"
    addEntrNHp = view(addEntr, nhmask)
    # ongoing part of RoME.jl #244
    addEntropyOnManifold!(mani, addEntrNHp, 1:getDimension(mani), spreadDist, pvec)
    # do info per coords
    ipc_ = zeros(getDimension(Xi[sfidx]))
    ipc_[pvec] .= 1.0
    ipc_
  end

  # check partial is easy as this is a prior
  return addEntr, ipc
end


function evalPotentialSpecific( Xi::AbstractVector{<:DFGVariable},
                                ccwl::CommonConvWrapper{Mixture{N_,F,S,T}},
                                solvefor::Symbol,
                                measurement::Tuple=(Vector{Vector{Float64}}(),);
                                kw... ) where {N_,F<:AbstractFactor,S,T}
  #
  evalPotentialSpecific(Xi,
                        ccwl,
                        solvefor,
                        F,
                        measurement;
                        kw... )
end

function evalPotentialSpecific( Xi::AbstractVector{<:DFGVariable},
                                ccwl::CommonConvWrapper{F},
                                solvefor::Symbol,
                                measurement::Tuple=(Vector{Vector{Float64}}(),);
                                kw... ) where {F <: AbstractFactor}
  #
  evalPotentialSpecific(Xi,
                        ccwl,
                        solvefor,
                        F,
                        measurement;
                        kw... )
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
                    skipSolve::Bool=false,
                    _slack=nothing  )
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
                                inflateCycles=inflateCycles, skipSolve=skipSolve,
                                _slack=_slack  )
  #
end

"""
    $SIGNATURES

Perform factor evaluation to resolve the "solve for" variable of a factor.  
This temporary function can be run without passing a factor graph object, but will internally allocate a new temporary new one.
Alternatively, the factor graph used for calculations can be passed in via the keyword `tfg`, hence the function name bang.

Notes
- `TypeParams_args::Vector{Tuple{InferenceVariable, P}}
- idea is please find best e.g. `b`, given `f(z,a,b,c)` either by roots or minimize (depends on definition of `f`)
- `sfidx::Int` is the solve for index, assuming `getVariableOrder(fct)`.

Example
```julia
B = _evalFactorTemporary!(EuclidDistance, (ContinuousScalar, ContinuousScalar), 2, ([10;],), ([0.],[9.5]) )
# should return `B = 10`
```

Related

[`evalFactor`](@ref), [`calcFactorResidual`](@ref), [`testFactorResidualBinary`](@ref)
"""
function _evalFactorTemporary!( fct::AbstractFactor,
                                varTypes::Tuple,
                                sfidx::Int,  # solve for index, assuming variable order for fct
                                meas_single::Tuple,
                                pts::Tuple;
                                tfg::AbstractDFG=initfg(),
                                solveKey::Symbol=:default,
                                newFactor::Bool=true,
                                _slack=nothing,
                                buildgraphkw... )
  #

  # build up a temporary graph in dfg
  _, _dfgfct = IIF._buildGraphByFactorAndTypes!(fct, varTypes, pts; dfg=tfg, solveKey=solveKey, newFactor=newFactor, buildgraphkw...)
  
  # get label convention for which variable to solve for 
  solvefor = getVariableOrder(_dfgfct)[sfidx]

  # do the factor evaluation
  measurement = ([meas_single[1],],)
  sfPts, _ = evalFactor(tfg, _dfgfct, solvefor, measurement, needFreshMeasurements=false, solveKey=solveKey, inflateCycles=1, _slack=_slack )

  # @info "EVALTEMP" length(sfPts) string(sfPts) meas_single measurement solvefor

  return sfPts
end