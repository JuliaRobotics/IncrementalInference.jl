
export calcFactorResidual



function approxConvBelief(dfg::AbstractDFG,
                          fc::DFGFactor,
                          target::Symbol,
                          measurement::AbstractVector=Tuple[];
                          solveKey::Symbol=:default,
                          N::Int=length(measurement), 
                          skipSolve::Bool=false )
  #
  v1 = getVariable(dfg, target)
  N = N == 0 ? getNumPts(v1, solveKey=solveKey) : N
  # points and infoPerCoord
  pts, ipc = evalFactor(dfg, fc, v1.label, measurement, solveKey=solveKey, N=N, skipSolve=skipSolve)
  
  len = length(ipc)
  mask = 1e-14 .< abs.(ipc)
  partl = collect(1:len)[ mask ]

  # is the convolution infoPerCoord full or partial
  if sum(mask) == len
    # not partial
    return manikde!(getManifold(getVariable(dfg, target)), pts, partial=nothing)
  else
    # is partial
    return manikde!(getManifold(getVariable(dfg, target)), pts, partial=partl)
  end
end


approxConv(w...;kw...) = getPoints( approxConvBelief(w...;kw...), false)

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
function approxConvBelief(dfg::AbstractDFG, 
                          from::Symbol, 
                          target::Symbol,
                          measurement::AbstractVector=Tuple[];
                          solveKey::Symbol=:default,
                          N::Int = length(measurement),
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
    @assert path[1] == from "sanity check failing for shortest path function"

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
    getBelief(dfg, varLbls[1]) |> getPoints
  else
    # chain would start one later
    idxS += 1
    # get the factor
    fct0 = getFactor(dfg,from)
    # get the Matrix{<:Real} of projected points
    pts1Bel = approxConvBelief(dfg, fct0, path[2], measurement, solveKey=solveKey, N=N, skipSolve=skipSolve)
    if length(path) == 2
      return pts1Bel
    end
    getPoints(pts1Bel)
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
      ptsBel = approxConvBelief(tfg, fct, path[idx+1], solveKey=solveKey, N=N, skipSolve=skipSolve)
      initManual!(tfg, path[idx+1], ptsBel)
      !setPPE ? nothing : setPPE!(tfg, path[idx+1], solveKey, ppemethod)
    end
  end

  # return target variable values
  return getBelief(tfg, target)
end



"""
    $(SIGNATURES)

Compute proposal belief on `vertid` through `fct` representing some constraint in factor graph.
Always full dimension variable node -- partial constraints will only influence subset of variable dimensions.
The remaining dimensions will keep pre-existing variable values.

Notes
- fulldim is true when "rank-deficient" -- TODO swap to false (or even float)
"""
function calcProposalBelief(dfg::AbstractDFG,
                            fct::DFGFactor,
                            target::Symbol,
                            measurement::AbstractVector=Tuple[];
                            N::Int=length(measurement),
                            solveKey::Symbol=:default,
                            dbg::Bool=false  )
  #
  # assuming it is properly initialized TODO
  proposal = approxConvBelief(dfg, fct, target, measurement, solveKey=solveKey, N=N)

  # return the proposal belief and inferdim, NOTE likely to be changed
  return proposal
end


# specifically the PartialPriorPassThrough dispatch
function calcProposalBelief(dfg::AbstractDFG,
                            fct::DFGFactor{<:CommonConvWrapper{<:PartialPriorPassThrough}},
                            target::Symbol,
                            measurement::AbstractVector=Tuple[];
                            N::Int=length(measurement),
                            solveKey::Symbol=:default,
                            dbg::Bool=false  )
  #

  # density passed through directly from PartialPriorPassThrough.Z
  fctFnc = getFactorType(fct)
  proposal = fctFnc.Z.densityFnc

  # in case of partial, place the proposal into larger marginal/partial MKD
  proposal_ = if isPartial(fctFnc)
    # oldbel = getBelief(dfg, target, solveKey)
    varType = getVariableType(dfg, target)
    M = getManifold(varType)
    u0 = getPointIdentity(varType)
    # replace(oldbel, proposal)
    antimarginal(M,u0,proposal,Int[fctFnc.partial...])
  else
    proposal
  end

  # return the proposal belief and inferdim, NOTE likely to be changed
  return proposal_
end

"""
    $SIGNATURES

Compute the proposals of a destination vertex for each of `factors` and place the result
as belief estimates in both `dens` and `partials` respectively.

Notes
- TODO: also return if proposals were "dimension-deficient" (aka ~rank-deficient).
"""
function proposalbeliefs!(dfg::AbstractDFG,
                          destlbl::Symbol,
                          factors::AbstractVector{<:DFGFactor},
                          dens::AbstractVector{<:ManifoldKernelDensity},
                          measurement::AbstractVector=Tuple[];
                          solveKey::Symbol=:default,
                          N::Int=getSolverParams(dfg).N, #maximum([length(getPoints(getBelief(dfg, destlbl, solveKey))); getSolverParams(dfg).N]),
                          dbg::Bool=false  )
  #

  # populate the full and partial dim containers
  ipcs = Vector{Vector{Float64}}(undef, length(factors))
  
  vardim = getDimension(getVariable(dfg, destlbl))
  # get a proposal belief from each factor connected to destlbl
  for (count,fct) in enumerate(factors)
    ccwl = _getCCW(fct)
    # need way to convey partial information
    # determine if evaluation is "dimension-deficient" solvable dimension
    # FIXME, update to infoPerCoord
    fct_ipc = ones(vardim) # getFactorSolvableDim(dfg, fct, destlbl, solveKey)
    # convolve or passthrough to get a new proposal
    propBel_ = calcProposalBelief(dfg, fct, destlbl, measurement, N=N, dbg=dbg, solveKey=solveKey)
    # partial density
    propBel = if isPartial(ccwl)
      pardims = _getDimensionsPartial(ccwl)
      @assert [getFactorType(fct).partial...] == [pardims...] "partial dims error $(getFactorType(fct).partial) vs $pardims"
      AMP.marginal(propBel_, Int[pardims...])
    else
      propBel_
    end
    push!(dens, propBel)
    ipcs[count] = fct_ipc
  end
  # len = maximum(length.(ipcs))
  ipc = zeros(vardim)
  for _ipc in ipcs
    ipc .+= _ipc
  end

  ipc
end
# group partial dimension factors by selected dimensions -- i.e. [(1,)], [(1,2),(1,2)], [(2,);(2;)]


# WIP, see `_buildGraphByFactorAndTypes!` where pts are full MKD Beliefs, following #1351 
# Legacy use in RoMEPlotting: plotFactor
# function approxConvBelief(fct::AbstractFactorRelative,
#                           varTypes::Union{<:Tuple,<:AbstractVector{<:InstanceType{T}}}, 
#                           mkds::Union{<:Tuple,<:AbstractVector{<:InstanceType{T}}};
#                           tfg::AbstractDFG=_buildGraphByFactorAndTypes!(fct,)
#                           ) where {T <: InferenceVariable}
#   #

# end



#
