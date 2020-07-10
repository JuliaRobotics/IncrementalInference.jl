


##==============================================================================
## Delete at end v0.13.x
##==============================================================================

export putMsgUpInitStatus!
export setTreeCliquesMarginalized!
export upPrepOutMsg!
export getCliqStatusUp, getCliqueStatusUp
export setCliqStatus!, getCliqStatus
export getMsgsUpChildrenInitDict
export doCliqUpSolve!
export fetchAssignTaskHistoryAll!, fetchCliqTaskHistoryAll!
export setCliqUpInitMsgs!
export getCliqInitUpMsgs, getInitDownMsg
export setMsgUpThis!, getMsgsUpThis
export setMsgDwnThis!, getMsgsDwnThis


function putMsgUpInitStatus!(cliq::TreeClique, status::CliqStatus, logger=SimpleLogger(stdout))
  @warn "putMsgUpInitStatus! is deprecated, try using prepPutCliqMsgUp! instead."
  cdat = getCliqueData(cliq)
  cdc = getMsgUpInitChannel_(cdat)
  cond = getSolveCondition(cliq)
    if isready(cdc)
      content = take!(cdc)
    end
  # FIXME, lock should not be required in all cases.
  lockUpStatus!(cliq, cliq.index, true, logger, true, "putMsgUpInitStatus!")
  setCliqueStatus!(cdat, status)
  put!(cdc, LikelihoodMessage(status=status))
  notify(cond)
    # FIXME hack to avoid a race condition  -- remove with atomic lock logic upgrade
    sleep(0.1)
    notify(cond) # getSolveCondition(cliq)
  #
  unlockUpStatus!(cdat)
  nothing
end

"""
    $SIGNATURES

Set all Bayes (Junction) tree cliques that have all marginalized and initialized variables.
"""
function setTreeCliquesMarginalized!(dfg::AbstractDFG,
                                     tree::AbstractBayesTree,
                                     logger=SimpleLogger(stdout))
  #
  for (cliid, cliq) in getCliques(tree)
    if areCliqVariablesAllMarginalized(dfg, cliq)

      ## FIXME, change to prepPutCliqueStatusMsgUp!
        # need to set the upward messages
        msgs = prepCliqInitMsgsUp(dfg, cliq)
        putMsgUpThis!(cliq, msgs)
        # TODO must be converted to pull model #674
        # prnt = getParent(tree, cliq)
        # if length(prnt) > 0
          # THIS IS FOR INIT PASSES ONLY
          putMsgUpInit!(cliq, msgs, logger)
        # end
        setCliqueStatus!(cliq, :marginalized)

      # set marginalized color
      setCliqDrawColor(cliq, "blue")

      # set flag, looks to be previously unused???
      getCliqueData(cliq).allmarginalized = true
    end
  end
  nothing
end

"""
    $SIGNATURES

Update clique status and notify of the change

Notes
- Assumes users will lock the status state before getting status until after decision whether to update status.
- If so, only unlock after status and condition has been updated.

Dev Notes
- Should be made an atomic transaction
"""
function notifyCliqUpInitStatus!(cliq::TreeClique,
                                 status::Symbol;
                                 logger=ConsoleLogger() )
  #
  @warn "notifyCliqUpInitStatus! is deprecated, use putMsgUpInitStatus! directly."
  cd = getCliqueData(cliq)
  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt), cliq=$(cliq.index), notifyCliqUpInitStatus! -- pre-lock, $(cd.initialized)-->$(status)"
  end
  flush(logger.stream)

  # currently using a lock internally (hack message channels are consolidated)
  putMsgUpInitStatus!(cliq, status, logger)

  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$(tt), cliq=$(cliq.index), notifyCliqUpInitStatus! -- unlocked, $(cd.initialized)"
  end

  nothing
end


"""
$(TYPEDEF)
"""
mutable struct MsgPassType
  fg::GraphsDFG
  cliq::TreeClique
  vid::Symbol # Int
  msgs::Array{LikelihoodMessage,1}
  N::Int
end


"""
    $SIGNATURES

Consolidation likely

DevNotes
- consolidation likely (prepCliqInitMsgsUp)
"""
function upPrepOutMsg!(dict::Dict{Symbol,TreeBelief}, seps::Vector{Symbol}, status::Symbol=:NULL)
  @error "upPrepOutMsg! is deprecated, use prepCliqInitMsgUP! instead."
  msg = LikelihoodMessage(status)
  for vid in seps
    msg.belief[vid] = dict[vid]
  end
  return msg
end

@deprecate getCliqueStatusUp(x...) getCliqueStatus(x...)
@deprecate getCliqStatusUp(x...) getCliqueStatus(x...)
@deprecate getCliqStatus(x...) getCliqueStatus(x...)
@deprecate setCliqStatus!(x...) setCliqueStatus!(x...)

@deprecate getMsgsUpChildrenInitDict(treel::AbstractBayesTree,cliq::TreeClique,::Type{TreeBelief},skip::Vector{Int}=Int[]) getMsgsUpInitChildren(treel, cliq, TreeBelief, skip)

@deprecate getMsgsUpChildrenInitDict(csmc::CliqStateMachineContainer,::Type{TreeBelief}=TreeBelief,skip::Vector{Int}=Int[] ) getMsgsUpInitChildren(csmc,TreeBelief,skip=skip )

@deprecate putMsgUpInit!(cliq::TreeClique,childid::Int,msg::LikelihoodMessage,logger=SimpleLogger(stdout)) putMsgUpInit!(cliq,msg,logger=logger)

@deprecate setMsgUpThisInitDict!(cdat::BayesTreeNodeData, idx, msg::LikelihoodMessage) setMsgUpThisInit!(cdat, msg)


"""
$(TYPEDEF)
"""
mutable struct UpReturnBPType
  upMsgs::LikelihoodMessage
  dbgUp::DebugCliqMCMC
  IDvals::Dict{Symbol, TreeBelief}
  keepupmsgs::LikelihoodMessage # Dict{Symbol, BallTreeDensity} # TODO Why separate upMsgs? FIXME consolidate with upMsgs
  totalsolve::Bool
  UpReturnBPType() = new()
  UpReturnBPType(x1,x2,x3,x4,x5) = new(x1,x2,x3,x4,x5)
end

@deprecate updateFGBT!(fg::AbstractDFG,cliq::TreeClique,urt::UpReturnBPType;dbg::Bool=false,fillcolor::String="",logger=ConsoleLogger()  ) updateFGBT!(fg, cliq, urt.IDvals, dbg=dbg, fillcolor=fillcolor, logger=logger)

@deprecate updateFGBT!(fg::AbstractDFG,bt::AbstractBayesTree,cliqID::Int,urt::UpReturnBPType;dbg::Bool=false, fillcolor::String=""  ) updateFGBT!( fg, getClique(bt, cliqID), urt, dbg=dbg, fillcolor=fillcolor )


function approxCliqMarginalUp!(fg_::AbstractDFG,
                               tree_::AbstractBayesTree,
                               cliq::TreeClique,
                               childmsgs=getMsgsUpChildren(fg_, tree_, cliq, TreeBelief);
                               N::Int=100,
                               dbg::Bool=false,
                               iters::Int=3,
                               drawpdf::Bool=false,
                               multiproc::Bool=true,
                               logger=ConsoleLogger()  )
  #
  error("OBSOLETE: approxCliqMarginalUp!(::AbstractDFG,...)")
  # approxCliqMarginalUp!(csmc.cliqSubFg, csmc.tree, csmc.cliq, getMsgsUpChildren(csmc, TreeBelief),N=N, dbg=dbg, iters=iters, drawpdf=drawpdf, multiproc=multiproc, logger=logger)
end

"""
    $SIGNATURES

Update `subfg<:AbstractDFG` according to internal computations for a full upsolve.
"""
function doCliqUpSolve!(csmc::CliqStateMachineContainer;
                        multiproc::Bool=getSolverParams(csmc.cliqSubFg).multiproc,
                        logger=ConsoleLogger()  )
  #
  approxCliqMarginalUp!(csmc, multiproc=multiproc, logger=logger)
  # csym = getCliqFrontalVarIds(csmc.cliq)[1]
  # approxCliqMarginalUp!(csmc, csym, false, N=getSolverParams(csmc.cliqSubFg).N, logger=logger, multiproc=multiproc)

  # TODO replace with msg channels only
  getCliqueData(csmc.cliq).upsolved = true
  return :upsolved
end

function approxCliqMarginalUp!(fgl::AbstractDFG,
                               treel::AbstractBayesTree,
                               csym::Symbol,
                               onduplicate::Bool;  # this must be deprecated for simplicity!
                               N::Int=100,
                               dbg::Bool=false,
                               iters::Int=3,
                               drawpdf::Bool=false,
                               multiproc::Bool=true,
                               logger=ConsoleLogger()  )
  #
  @warn "approxCliqMarginalUp! API is changing, use csmc version instead."
  @assert !onduplicate "approxCliqMarginalUp! onduplicate keyword is being deprecated"
  fg_ = onduplicate ? deepcopy(fgl) : fgl
  # onduplicate
  with_logger(logger) do
    @warn "rebuilding new Bayes tree on deepcopy of factor graph"
  end
  # FIXME, should not be building a new tree here since variable orderings can differ!!!
  tree_ = onduplicate ? wipeBuildNewTree!(fgl) : treel

  # copy up and down msgs that may already exists #TODO Exists where? it copies from tree_ to tree
  if onduplicate
    for (id, cliq) in treel.cliques
      setUpMsg!(tree_.cliques[cliq.index], getUpMsgs(cliq)) #TODO cliq.index may be problematic, how do we know it will be the same index on rebuilding?
      setDwnMsg!(tree_.cliques[cliq.index], getDwnMsgs(cliq))
    end
  end

  cliq = getCliq(tree_, csym)
  # setCliqDrawColor(cliq, "red")

  approxCliqMarginalUp!(fg_, tree_, cliq, N=N, dbg=dbg, iters=iters, drawpdf=drawpdf, multiproc=multiproc, logger=logger)
end

function doCliqUpSolve!(subfg::AbstractDFG,
                        tree::AbstractBayesTree,
                        cliq::TreeClique;
                        multiproc::Bool=true,
                        logger=ConsoleLogger()  )
  #
  @error "doCliqUpSolve! is being refactored, use the csmc version instead"
  csym = getCliqFrontalVarIds(cliq)[1]
  # csym = DFG.getVariable(subfg, getCliqFrontalVarIds(cliq)[1]).label # ??
  approxCliqMarginalUp!(subfg, tree, csym, false, N=getSolverParams(subfg).N, logger=logger, multiproc=multiproc)
  getCliqueData(cliq).upsolved = true
  return :upsolved
end

function doCliqInitDown!(subfg::AbstractDFG,
                         tree::AbstractBayesTree,
                         cliq::TreeClique;
                         dbg::Bool=false )
  #
  @error("deprecated doCliqInitDown!(subfg, tree, cliq) use doCliqInitDown!(subfg, cliq, dwinmsgs) instead.")
  prnt = getParent(tree, cliq)[1]
  dwinmsgs = prepCliqInitMsgsDown!(subfg, tree, prnt)
  status = doCliqInitDown!(subfg, cliq, dwinmsgs, dbg=dbg)

  return status
end

@deprecate fetchCliqTaskHistoryAll!(x...) fetchCliqHistoryAll!(x...)

function fetchAssignTaskHistoryAll!(tree::AbstractBayesTree, smt)
  hist = Dict{Int, Vector{Tuple{DateTime,Int,Function,CliqStateMachineContainer}}}()
  fetchCliqTaskHistoryAll!(smt, hist)
  assignTreeHistory!(tree, hist)
end

@deprecate getMsgsUpChildren(::AbstractDFG, treel::AbstractBayesTree, cliq::TreeClique, ::Type{TreeBelief}) getMsgsUpChildren(treel,cliq,TreeBelief)


@deprecate setCliqUpInitMsgs!(x...) putMsgUpInit!(x...)

@deprecate getCliqInitUpMsgs(x...) getMsgUpThisInit(x...)
@deprecate getInitDownMsg(x...) getMsgDwnThisInit(x...)

@deprecate getMsgsUpThis(x...) fetchMsgUpThis(x...)
@deprecate setMsgUpThis!(x...) putMsgUpThis!(x...)

@deprecate getMsgsDwnThis(x...) fetchMsgDwnThis(x...)
@deprecate setMsgDwnThis!(x...) putMsgDwnThis!(x...)



##==============================================================================
## Delete at end v0.12.x
##==============================================================================

# export getCliqparentMsgDown
export setDwnMsg!
export upMsg, dwnMsg
export getDwnMsgs
export getCliq, whichCliq, hasCliq
export getCliqChildMsgsUp
export setUpMsg!, getUpMsgs
export assignTreeHistory!
export getVertKDE,  getVert


# # return ::Vector{DFGFactor}
# # TODO, perhaps consolidate
# function addMsgFactors_Parametric!(subfg::AbstractDFG,
#                                    msgs::LikelihoodMessage)
#   # add messages as priors to this sub factor graph
#   msgfcts = DFGFactor[]
#   svars = DFG.listVariables(subfg)
#   for (msym, belief_) in msgs.belief
#     if msym in svars
#       #TODO covaraince
#       #TODO Maybe always use MvNormal
#       if size(belief_.val, 2) == 1 && size(belief_.val, 1) == 1
#         msgPrior =  MsgPrior(Normal(belief_.val[1], sqrt(belief_.bw[1])), belief_.inferdim)
#       elseif size(belief_.val, 2) == 1 && 1 < size(belief_.val, 1)
#         mvnorm = createMvNormal(belief_.val[:,1], belief_.bw)
#         mvnorm != nothing ? nothing : (return DFGFactor[])
#         msgPrior =  MsgPrior(mvnorm, belief_.inferdim)
#       else
#         error("Don't know what what to do with size(belief_.val)=$(size(belief_.val))")
#       end
#       fc = addFactor!(subfg, [msym], msgPrior, graphinit=false)
#       push!(msgfcts, fc)
#     end
#   end
#   return msgfcts
# end

# # consolidated with AbstractBayesTree in TreeMessageAccessors.jl
# function takeBeliefMessageUp!(tree::MetaBayesTree, edge)
#   # Blocks until data is available.
#   beliefMsg = take!(getMsgUpChannel(tree, edge))
#   return beliefMsg
# end
# function takeBeliefMessageDown!(tree::MetaBayesTree, edge)
#   # Blocks until data is available.
#   beliefMsg = take!(getMsgDwnChannel(tree, edge))
#   return beliefMsg
# end
# function putBeliefMessageDown!(tree::MetaBayesTree, edge, beliefMsg::LikelihoodMessage)
#   # Blocks until data is available.
#   put!(getMsgDwnChannel(tree, edge), beliefMsg)
#   return beliefMsg
# end
# function putBeliefMessageUp!(tree::MetaBayesTree, edge, beliefMsg::LikelihoodMessage)
#   # Blocks until data is available.
#   put!(getMsgUpChannel(tree, edge), beliefMsg)
#   return beliefMsg
# end


# better version in TreeMessageUtils.jl
# function getMsgsUpChildren(fg_::AbstractDFG,
#                            treel::AbstractBayesTree,
#                            cliq::TreeClique,
#                            ::Type{TreeBelief} )
  #
  # childmsgs = LikelihoodMessage[]
  # for child in getChildren(treel, cliq)
  #   nbpchild = LikelihoodMessage()
  #   for (key, bel) in getUpMsgs(child).belief
  #     # manis = getManifolds(fg_, key)
  #     # inferdim = getVariableInferredDim(fg_, key)
  #     dcBel = deepcopy(bel)
  #     nbpchild.belief[key] = TreeBelief(dcBel.val, dcBel.bw, dcBel.inferdim, getSofttype(getVariable(fg_, key)))
  #   end
  #   push!(childmsgs, nbpchild)
  # end
  # return childmsgs
# end

function getMsgsUpChildren(treel::AbstractBayesTree,
                            cliq::TreeClique,
                            ::Type{BallTreeDensity})
  #
  @warn "BallTreeDensity version of getMsgsUpChildren is deprecated, use TreeBelief version instead."
  childmsgs = IntermediateMultiSiblingMessages()
  for child in getChildren(treel, cliq)
    for (key, bel) in getUpMsgs(child).belief
      # id = fg_.IDs[key]
      # manis = getManifolds(fg_, id)
      if !haskey(childmsgs, key)
        childmsgs[key] = IntermediateSiblingMessages()
      end
      push!(childmsgs[key], bel )
    end
  end
  return childmsgs
end

# TODO consolidate
function getMsgsUpChildren(csmc::CliqStateMachineContainer,
                            ::Type{BallTreeDensity})
  #
  @warn "BallTreeDensity version of getMsgsUpChildren is deprecated, use TreeBelief version instead."
  getMsgsUpChildren(csmc.tree, csmc.cliq, BallTreeDensity)
end

function addMsgFactors!(subfg::AbstractDFG,
                        msgs::Dict{Symbol, Vector{Tuple{BallTreeDensity, Float64}}} )
      # msgs::
  # add messages as priors to this sub factor graph
  @warn "Tuple{KDE,Floa64} specific version of addMsgFactors! is deprecated, use LikelihoodMessage version instead."
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, dms) in msgs
    for dm in dms
      if msym in svars
        # TODO should be on manifold prior, not just generic euclidean prior -- okay since variable on manifold, but not for long term
        fc = addFactor!(subfg, [msym], MsgPrior(dm[1], dm[2]), graphinit=false)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end


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


#
