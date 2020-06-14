


##==============================================================================
## Delete at end v0.12.x
##==============================================================================

export getCliqparentMsgDown
export setDwnMsg!
export upMsg, dwnMsg
export getDwnMsgs
export getCliq, whichCliq, hasCliq
export getCliqChildMsgsUp
export setUpMsg!, getUpMsgs
export assignTreeHistory!
export getVertKDE,  getVert

@deprecate LikelihoodMessage(status::CliqStatus) LikelihoodMessage(status=status)
@deprecate LikelihoodMessage(status::CliqStatus, varOrder::Vector{Symbol}, cliqueLikelihood::SamplableBelief) LikelihoodMessage(status=status, variableOrder=varOrder, cliqueLikelihood=cliqueLikelihood)
@deprecate LikelihoodMessage(status::CliqStatus, cliqueLikelihood::SamplableBelief) LikelihoodMessage(status=status, cliqueLikelihood=cliqueLikelihood)

"""
    $SIGNATURES

Build a new subgraph from `fgl<:AbstractDFG` containing all variables and factors
associated with `cliq`.  Additionally add the upward message prior factors as
needed for belief propagation (inference).

Notes
- `cliqsym::Symbol` defines the cliq where variable appears as a frontal variable.
- `varsym::Symbol` defaults to the cliq frontal variable definition but can in case a
  separator variable is required instead.
"""
function buildCliqSubgraphDown(fgl::AbstractDFG, treel::AbstractBayesTree, cliqsym::Symbol, varsym::Symbol=cliqsym)
  @warn "Obsolete, buildCliqSubGraph*() is no longer in use"
  # build a subgraph copy of clique
  cliq = whichCliq(treel, cliqsym)
  syms = getCliqAllVarIds(cliq)
  subfg = buildSubgraph(fgl, syms, 1)

  # add upward messages to subgraph
  msgs = getMsgDownParent(treel, cliq)
  addMsgFactors!(subfg, msgs)
  return subfg
end


@deprecate getCliqParentMsgDown(x...) getMsgDwnParent(x...)

# getCliq(bt::AbstractBayesTree, frt::Symbol) = getClique(bt, bt.frontals[frt])
# whichCliq(bt::AbstractBayesTree, frt::Symbol) = getCliq(bt, frt)
# whichCliq(bt::AbstractBayesTree, frt::AbstractString) = whichCliq(bt, Symbol(frt))

@deprecate getCliq(x...) getClique(x...)
@deprecate whichCliq(x...) getClique(x...)
@deprecate hasCliq(x...) hasClique(x...)

@deprecate getCliqChildMsgsUp(x...) getMsgsUpChildren(x...)

# export getCliqPotentials
# @deprecate getCliqPotentials(dfg::AbstractDFG,bt::AbstractBayesTree,cliq::TreeClique) getCliquePotentials(dfg, bt, cliq)

@deprecate upMsg(x...) getMsgsUpThis(x...)
@deprecate dwnMsg(x...) getMsgsDwnThis(x...)
@deprecate getDwnMsgs(x...) getMsgsDwnThis(x...)
@deprecate setDwnMsg!(x...) setMsgDwnThis!(x...)
@deprecate setUpMsg!(cliql::TreeClique, msgs::LikelihoodMessage) setMsgUpThis!(cliql, msgs)
@deprecate getUpMsgs(x...) getMsgsUpThis(x...)

# NOTE decided not to store messages in CSMC, but closer to Tree instead (likely on edges)
# function setUpMsg!(csmc::CliqStateMachineContainer, cliqid::Int, msgs::LikelihoodMessage)
#   csmc.msgsUp[cliqid] = msgs
# end
# getUpMsgs(csmc::CliqStateMachineContainer) = csmc.msgsUp

"""
    $SIGNATURES

Return clique state machine history from `tree` if it was solved with `recordcliqs`.

Notes
- Cliques are identified by front variable `::Symbol` which are always unique across the cliques.
"""
function getCliqSolveHistory(cliq::TreeClique)
  @error ".statehistory is obsolete, use fetch.(smt) instead."
  # getCliqueData(cliq).statehistory
end
function getCliqSolveHistory(tree::AbstractBayesTree, frntal::Symbol)
  cliq = whichCliq(tree, frntal)
  getCliqSolveHistory(cliq)
end

"""
    $SIGNATURES

Return dict of all histories in a Bayes Tree.
"""
function getTreeCliqsSolverHistories(fg::G,
                                     tree::AbstractBayesTree)::Dict{Symbol, CSMHistory} where G <: AbstractDFG
  #
  @error "obsolete"
  # fsy = getTreeAllFrontalSyms(fg, tree)
  # histories = Dict{Symbol, CSMHistory}()
  # for fs in fsy
  #   hist = getCliqSolveHistory(tree, fs)
  #   if length(hist) > 0
  #     histories[fs] = hist
  #   end
  # end
  # return histories
end

function printCliqHistorySummary(cliq::TreeClique)
  hist = getCliqSolveHistory(cliq)
  printCliqHistorySummary(hist)
end

function printCliqHistorySummary(tree::AbstractBayesTree, frontal::Symbol)
  hist = getCliqSolveHistory(tree, frontal)
  printCliqHistorySummary(hist)
end

"""
    $SIGNATURES

After solving, clique histories can be inserted back into the tree for later reference.
This function helps do the required assigment task.
"""
function assignTreeHistory!(treel::AbstractBayesTree, cliqHistories::Dict)
  @error "assignTreeHistory! is obsolete."
  # for i in 1:length(getCliques(treel))
  #   if haskey(cliqHistories, i)
  #     hist = cliqHistories[i]
  #     for i in 1:length(hist)
  #       hist[i][4].logger = SimpleLogger(stdout)
  #     end
  #     # getCliqueData(treel, i).statehistory=hist
  #   end
  # end
end


@deprecate emptyBTNodeData() BayesTreeNodeData()


function evalPotentialSpecific(Xi::Vector{DFGVariable},
                               ccwl::CommonConvWrapper{T},
                               solvefor::Symbol,
                               measurement::Tuple=(zeros(0,100),);
                               N::Int=size(measurement[1],2),
                               spreadNH::Real=3.0,
                               dbg::Bool=false ) where {T <: FunctorPairwiseNH}
  #
  @warn "FunctorPairwiseNH will be deprecated in favor of common `nullhypo=` interface."
  # TODO -- could be constructed and maintained at addFactor! time
  sfidx, maxlen, manis = prepareCommonConvWrapper!(ccwl, Xi, solvefor, N)
  # prepare nullhypothesis
  allelements, nhc, ENT = assembleNullHypothesis(ccwl, maxlen, spreadNH)

  # Compute across the true or null hypothesis
  computeAcrossNullHypothesis!(ccwl, allelements, nhc, ENT )

  return ccwl.params[ccwl.varidx]
end


function evalPotentialSpecific(Xi::Vector{DFGVariable},
                               ccwl::CommonConvWrapper{T},
                               solvefor::Symbol,
                               measurement::Tuple=(zeros(0,100),);
                               N::Int=size(measurement[1],2),
                               spreadfactor::Real=10.0,
                               dbg::Bool=false,
                               spreadNH::Float64=3.0 ) where {T <: FunctorSingletonNH}
  #
  @warn "FunctorSingletonNH will be deprecated in favor of common `nullhypo=` interface."
  fnc = ccwl.usrfnc!

  val = getVal(Xi[1])
  d = size(val,1)
  var = Statistics.var(val, dims=2) .+ 1e-3

  # prep in case special samplers used
  # determine amount share of null hypothesis particles
  freshSamples!(ccwl, N, FactorMetadata(), Xi)
  # ccwl.measurement = getSample(ccwl.usrfnc!, N)
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


function convert(::Type{TreeBelief},
                 bel::Tuple{BallTreeDensity,Float64},
                 manifolds::T) where {T <: Tuple}
  @error "Dont use this convert(::Type{TreeBelief}, bel::Tuple{BallTreeDensity,Float64}, manifolds) since it must assume ContinuousScalar softtype!!!"
  TreeBelief(getPoints(bel[1]), getBW(bel[1])[:,1:1], bel[2], ContinuousScalar(), manifolds)
end

"""
    $(SIGNATURES)

Encode complicated function node type to related 'Packed<type>' format assuming a user supplied convert function .
"""
function convert2packedfunctionnode(fgl::G,
                                    fsym::Symbol ) where G <: AbstractDFG
  #
  @warn "convert2packedfunctionnode is obsolete and will be removed, see DFG serialization."
  # fid = fgl.fIDs[fsym]
  fnc = getfnctype(fgl, fsym)
  usrtyp = convert(PackedInferenceType, fnc)
  cfnd = convert(PackedFunctionNodeData{usrtyp}, getSolverData(getFactor(fgl, fsym)) )
  return cfnd, usrtyp
end



@deprecate getVertKDE(v::DFGVariable) getKDE(v)
@deprecate getVertKDE(dfg::AbstractDFG, lbl::Symbol) getKDE(dfg, lbl)


function edgelist2edgedict(edgelist::Array{Graphs.Edge{TreeClique},1})
  error("edgelist2edgedict is obsolete, use DFG methods instead.")
  edgedict = Dict{Int,Graphs.Edge{TreeClique}}()
  for edge in edgelist
    edgedict[edge.index] = edge
  end
  return edgedict
end


# TODO: Confirm this is supposed to be a variable?
function setVal!(v::DFGVariable, em::TreeBelief; solveKey::Symbol=:default)
    @warn "setVal! deprecated, use setValKDE! instead"
    setValKDE!(v, em, solveKey=solveKey)
end
function setVal!(v::DFGVariable, p::BallTreeDensity; solveKey::Symbol=:default)
    @warn "setVal! deprecated, use setValKDE! instead"
    setValKDE!(v, p, solveKey=solveKey)
end

@deprecate manualinit!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity)
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) initManual!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) initManual!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) false

export _evalType

# not sure if and where this is still being used
function _evalType(pt::String)::Type
    @error "_evalType has been deprecated, use DFG serialization methods instead."
    try
        getfield(Main, Symbol(pt))
    catch ex
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
    end
end

"""
    $SIGNATURES

writeGraphPdf deprecated, use drawGraph instead
"""
function writeGraphPdf(fgl::AbstractDFG;
                       viewerapp::AbstractString="evince",
                       filepath::AbstractString="/tmp/fg.pdf",
                       engine::AbstractString="neato",
                       show::Bool=true )
  #
  @warn "writeGraphPdf deprecated, use drawGraph instead"
  drawGraph(fgl, viewerapp=viewerapp, filepath=filepath, engine=engine, show=show )
end

function getVert(dfg::AbstractDFG, sym::Symbol, nt::Symbol=:var)
  @warn "IIF.getVert is deprecated, use DFG.getVariable or DFG.getFactor instead."
  if nt == :var
    return DFG.getVariable(dfg, sym)
  elseif nt == :fct
    return DFG.getFactor(dfg, sym)
  else
    error("unknown getVert request nt=$nt")
  end
end

##==============================================================================
## Delete at end v0.11.x
##==============================================================================

#
