
##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

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


##==============================================================================
## Delete when tree message consolidation is complete
##==============================================================================


function messages(btnd::BayesTreeNodeData)
  @warn("btnd.messages will be deprecated")
  btnd.messages
end

messages(clique::TreeClique) = getCliqueData(clique).messages


##==============================================================================
## Delete at end v0.16.x
##==============================================================================

@deprecate wipeBuildNewTree!(dfg::AbstractDFG; kwargs...) resetBuildTree!(dfg; kwargs...)

# getSample(s::MixtureRelative, N::Int=1) = (reshape.(rand.(s.Z, N),1,:)..., rand(s.C, N))

@deprecate (MixturePrior{T}(z::NTuple{N,<:SamplableBelief}, c::Union{<:Distributions.DiscreteNonParametric, NTuple{N,<:Real}, <:AbstractVector{<:Real}} ) where {T,N}) MixturePrior(z,c)

@deprecate LinearConditional(N::Int=1) LinearRelative{N}(LinearAlgebra.I)
# @deprecate LinearConditional(x::SamplableBelief) LinearRelative(x)
@deprecate LinearConditional(x...) LinearRelative(x...)

# function LinearConditional{N, T}(x...) where {N,T}
#   @warn("LinearConditional{N, T} is deprecated, use LinearRelative instead")
#   LinearRelative{N,T}(x...)
# end

@deprecate PackedLinearConditional(x...) PackedLinearRelative(x...)



## =============================================================================
## Clique condition locks
## =============================================================================

@deprecate lockUpStatus!(x...) ()->nothing

@deprecate unlockUpStatus!(cdat::BayesTreeNodeData) ()->nothing
unlockUpStatus!(cliq::TreeClique) = unlockUpStatus!(getCliqueData(cliq))

@deprecate lockDwnStatus!(cdat::BayesTreeNodeData, idx::Int=1; logger=ConsoleLogger()) ()->nothing

@deprecate unlockDwnStatus!(cdat::BayesTreeNodeData) ()->nothing




function iifdepwarn(msg, funcsym; maxlog=nothing)
  @logmsg(
      Base.CoreLogging.Warn,
      msg,
      _module=begin
          bt = backtrace()
          frame, caller = Base.firstcaller(bt, funcsym)
          # TODO: Is it reasonable to attribute callers without linfo to Core?
          caller.linfo isa Core.MethodInstance ? caller.linfo.def.module : Core
      end,
      _file=String(caller.file),
      _line=caller.line,
      _id=(frame,funcsym),
      _group=:iifdepwarn,
      caller=caller,
      short_stacktrace=stacktrace(bt)[7:9],
      maxlog=maxlog
  )
  nothing
end

function Base.getproperty(obj::BayesTreeNodeData, sym::Symbol)
  if sym == :dwnMsg
    # iifdepwarn("#459 get dwnMsg", :getproperty)
  elseif sym == :downInitMsg
    # iifdepwarn("#459 get downInitMsg", :getproperty)
  elseif sym == :initDownChannel
    # iifdepwarn("#459 get initDownChannel", :getproperty)
  end
  return getfield(obj, sym)
end

function Base.setproperty!(obj::BayesTreeNodeData, sym::Symbol, val)
  if sym == :dwnMsg
    # iifdepwarn("#459 set dwnMsg", :setproperty!)
  elseif sym == :downInitMsg
    # iifdepwarn("#459 set downInitMsg", :setproperty!)
  elseif sym == :initDownChannel
    # iifdepwarn("#459 set initDownChannel", :setproperty!)
  end
  return setfield!(obj, sym, convert(fieldtype(typeof(obj), sym), val))
end


## ============================================================================
## .initDownChannel, MUST BE REMOVED
## ============================================================================

## ============================================================================
## .downInitMsg, MUST BE REMOVED
## ============================================================================


@deprecate putMsgDwnInitChannel!(btnd::BayesTreeNodeData, msg::LikelihoodMessage) putDwnMsgConsolidated!(btnd, msg)
@deprecate getMsgDwnInitChannel_(btnd::BayesTreeNodeData) getDwnMsgConsolidated(btnd)
# getMsgDwnInitChannel_(btnd::BayesTreeNodeData) = btnd.initDownChannel

function getMsgDwnInitChannel_(cliq::TreeClique)
  @warn("getMsgDwnInitChannel_ is deprecated, use getDwnMsgConsolidated instead.")
  getMsgDwnInitChannel_(getCliqueData(cliq))
end
fetchMsgDwnInit(cliq::TreeClique) = fetch(getMsgDwnInitChannel_(cliq))


# FIXME OLD must be consolidated as part of 459
function putMsgDwnInitStatus!(cliq::TreeClique, status::CliqStatus, logger=ConsoleLogger(), msg=LikelihoodMessage(status=status))
  @warn("putMsgDwnInitStatus! is deprecated, use putDwnMsgConsolidated! instead")
  cdat = getCliqueData(cliq)
  cdc = getDwnMsgConsolidated(cdat)
  # cdc = getMsgDwnInitChannel_(cdat)
    if isready(cdc)
      content = take!(cdc)
      with_logger(logger) do
        @info "dumping stale cliq=$(cliq.index) status message $(content), replacing with $(status)"
      end
    end
  put!(cdc, msg)
  notify(getSolveCondition(cliq))
    # FIXME hack to mitigate old race condition
    sleep(0.1)
    notify(getSolveCondition(cliq))

  nothing
end

function getfetchCliqueInitMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing)
  @debug "getfetchCliqueInitMsgDown from=$(from)"
  return cdata.downInitMsg
end
# getMsgDwnThisInit(cliq::TreeClique) = getMsgDwnThisInit(getCliqueData(cliq)) # WHAT ???

function putCliqueInitMsgDown!(cdata::BayesTreeNodeData, initmsg::LikelihoodMessage)
  @warn("putCliqueInitMsgDown! is deprecated, use putDwnMsgConsolidated instead")
  cdata.downInitMsg = initmsg
  nothing
end


@deprecate fetchMsgDwnThis(cliq::TreeClique) fetchDwnMsgConsolidated(cliq)

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
fetchMsgDwnThis(csmc::CliqStateMachineContainer) = fetchMsgDwnThis(csmc.cliq)
fetchMsgDwnThis(btl::AbstractBayesTree, sym::Symbol) = fetchMsgDwnThis(getClique(btl, sym))
# function fetchMsgDwnThis(cliql::TreeClique)
#   @error("fetchMsgDwnThis is deprecated, use getDwnMsgConsolidated instead")
#  # fetchDwnMsgConsolidated(cliql)
  # getCliqueData(cliql).dwnMsg
# end


@deprecate putMsgDwnThis!(x...; kw...) putDwnMsgConsolidated!(x...; kw...)

"""
$(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""  
putMsgDwnThis!(csmc::CliqStateMachineContainer, msgs::LikelihoodMessage) = putMsgDwnThis!(csmc.cliq, msgs)  # NOTE, old, csmc.msgsDown = msgs
# function putMsgDwnThis!(cdata::BayesTreeNodeData, msg::LikelihoodMessage; from::Symbol=:nothing)
#   @error("putMsgDwnThis! is deprecated, use putDwnMsgConsolidated instead")
#   @debug "putMsgDwnThis! from=$(from)"
#   cdata.dwnMsg = msg
# end  
# function putMsgDwnThis!(cliql::TreeClique, msgs::LikelihoodMessage)
#   getCliqueData(cliql).dwnMsg = msgs
# end  
  




"""
    $SIGNATURES

Do down initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in downward direction.

Notes
- State machine function nr. 8a
- Includes initialization routines.
- TODO: Make multi-core

DevNotes
- FIXME major refactor of this function required.
- FIXME this function actually occur during the parent CSM, therefore not all pull model #674
"""
function collectDwnInitMsgFromParent_StateMachine(csmc::CliqStateMachineContainer)
  #
  # TODO consider exit early for root clique rather than avoiding this function
  infocsm(csmc, "8a, needs down message -- attempt down init")
  setCliqDrawColor(csmc.cliq, "gold")

  # initialize clique in downward direction
  # not if parent also needs downward init message
  # prnt = getParent(csmc.tree, csmc.cliq)[1]
  opt = getSolverParams(csmc.dfg)
  @assert !haskey(opt.devParams,:dontUseParentFactorsInitDown) "dbgnew is old school, 459 dwninit consolidation has removed the option for :dontUseParentFactorsInitDown"

  # take atomic lock OF PARENT ??? when waiting for downward information
  # lockUpStatus!(prnt, prnt.index, true, csmc.logger, true, "cliq$(csmc.cliq.index)") # TODO XY ????
  # infocsm(csmc, "8a, after up lock")

  # get down message from the parent
  # check if any msgs should be multiplied together for the same variable
  # get the current messages ~~stored in~~ [going to] the parent (pull model #674)
  # FIXME, post #459 calls?
  # this guy is getting any sibling up messages by calling on the parent
  prntmsgs::Dict{Int, LikelihoodMessage} = getMsgsUpInitChildren(csmc.tree, csmc.cliq, TreeBelief, skip=[csmc.cliq.index;])         
  
  # reference to default dict location
  dwinmsgs = getfetchCliqueInitMsgDown(csmc.cliq.data, from=:getMsgDwnThisInit) |> deepcopy  #JT 459 products = getMsgDwnThisInit(prnt)
  infocsm(csmc, "getMsgInitDwnParent -- msg ids::Int=$(collect(keys(prntmsgs)))")
  
  # stack all parent incoming upward messages into dict of vector msgs
  prntBelDictVec::Dict{Symbol, Vector{TreeBelief}} = convertLikelihoodToVector(prntmsgs, logger=csmc.logger)
  ## TODO use parent factors too
  # intersect with the asking clique's separator variables
  # this function populates `dwinmsgs` with the appropriate products described in `prntBelDictVec`
  # FIXME, should not be using full .dfg ???
  condenseDownMsgsProductPrntFactors!(csmc.dfg, dwinmsgs, prntBelDictVec, prnt, csmc.cliq, csmc.logger)

  # remove msgs that have no data
  rmlist = Symbol[]
  for (prsym,belmsg) in dwinmsgs.belief
    if belmsg.inferdim < 1e-10
      # no information so remove
      push!(rmlist, prsym)
    end
  end
  infocsm(csmc, "prepCliqInitMsgsDown! -- rmlist, no inferdim, keys=$(rmlist)")
  for pr in rmlist
    delete!(dwinmsgs.belief, pr)
  end

  infocsm(csmc, "prepCliqInitMsgsDown! -- product keys=$(collect(keys(dwinmsgs.belief)))")

  # now put the newly computed message in the appropriate container
  # FIXME THIS IS A PUSH MODEL, see #674 -- must make pull model first
  # FIXME must be consolidated as part of #459
  putCliqueInitMsgDown!(getCliqueData(csmc.cliq), dwinmsgs)
  # unlock
  # unlockUpStatus!(prnt) # TODO XY ????, maybe can remove after pull model #674?
  infocsm(csmc, "8a, attemptCliqInitD., unlocked")

  # go to 7b (maybe, and part of dwnMsg #459 WIP 9)
  return slowIfChildrenNotUpSolved_StateMachine
    # # go to 8j.
    # return dwnInitSiblingWaitOrder_StateMachine
end


# currently for internal use only
# initialize variables based on best current achievable ordering
# OBVIOUSLY a lot of refactoring and consolidation needed with cliqGibbs / approxCliqMarginalUp
function initSolveSubFg!( subfg::AbstractDFG,
                          logger=ConsoleLogger() )
  #
  error("initSolveSubFg! is deprecated, do not use.  Follow post #459 CSM")
  varorder = getSubFgPriorityInitOrder(subfg, logger)
  with_logger(logger) do
    @info "initSolveSubFg! -- varorder=$varorder"
  end
  cycleInitByVarOrder!(subfg, varorder, logger=logger)
  nothing
end




# Helper function for prepCliqInitMsgsDown!
# future, be used in a cached system with parent in one location only for all siblings
# this function rebuilds a local subgraph from dfg and performs the calculations of the parent here and does not wait on the CSM to perform anything.
# 4-stroke compute may render this whole function obsolete.
function condenseDownMsgsProductPrntFactors!( fgl::AbstractDFG,
                                              products::LikelihoodMessage,
                                              msgspervar::Dict{Symbol, <:AbstractVector},
                                              prnt::TreeClique,
                                              cliq::TreeClique,
                                              logger=ConsoleLogger() )
  #
  error("condenseDownMsgsProductPrntFactors! is deprecated, follow post #459 instead")

  # determine the variables of interest
  reqMsgIds = collect(keys(msgspervar))
  # unique frontals per cliq
  prntvars = intersect(getCliqSeparatorVarIds(cliq), getCliqAllVarIds(prnt))
  lvarids = union(prntvars, reqMsgIds)
  # determine allowable factors, if any (only from parent cliq)
  awfcts = getCliqFactorIdsAll(prnt)

  # build required subgraph for parent/sibling down msgs
  lsfg = buildSubgraph(fgl, lvarids, 1; verbose=false)

  tempfcts = lsf(lsfg)
  dellist = setdiff(awfcts, tempfcts)
  for delf in dellist
    # TODO -- double check this deletefactor method is leaving the right parent sharing factor graph behind
    if exists(lsfg, delf)
      deleteFactor!(lsfg,delf)
    end
  end

  # add message priors
  addMsgFactors!(lsfg, msgspervar) # , DownwardPass

  # perform initialization/inference
  # FIXME, better consolidate with CSM ....uhhh TODO
  initSolveSubFg!(lsfg, logger)

  # extract complete downward marginal msg priors
  for id in intersect(getCliqSeparatorVarIds(cliq), lvarids)
    vari = getVariable(lsfg, id)
    products.belief[id] = TreeBelief(vari)
  end

  nothing
end
# # QUICK DBG CODE
# with_logger(logger) do
#     @info "condenseDownMsgsProductPrntFactors! -- reqMsgIds=$(reqMsgIds),"
#     @info "condenseDownMsgsProductPrntFactors! -- vars=$(lvarids),"
#     @info "condenseDownMsgsProductPrntFactors! -- allow factors $awfcts"
# end
# with_logger(logger) do
#     @info "condenseDownMsgsProductPrntFactors! -- lsfg fcts=$(lsf(lsfg)),"
#     @info "condenseDownMsgsProductPrntFactors! -- excess factors $dellist"
# end
# vars = ls(lsfg)
# len = length(vars)
# tdims = Vector{Float64}(undef, len)
# isinit = Vector{Bool}(undef, len)
# for idx in 1:len
#   tdims[idx] = getVariableSolvableDim(lsfg, vars[idx])
#   isinit[idx] = isInitialized(lsfg, vars[idx])
# end
# with_logger(logger) do
#     @info "condenseDownMsgsProductPrntFactors! -- after cycle init: vars=$vars"
#     @info "condenseDownMsgsProductPrntFactors! -- after cycle init: tdims=$tdims"
#     @info "condenseDownMsgsProductPrntFactors! -- after cycle init: isinit=$isinit"
# end
# # QUICK DBG CODE




"""
    $SIGNATURES

Initialization requires down message passing of more specialized down init msgs.
This function performs any possible initialization of variables and retriggers
children cliques that have not yet initialized.

Notes:
- Assumed this function is only called after status from child clique up inits completed.
- Assumes cliq has parent.
  - will fetch message from parent
- Will perform down initialization if status == `:needdownmsg`.
- might be necessary to pass furhter down messges to child cliques that also `:needdownmsg`.
- Will not complete cliq solve unless all children are `:upsolved` (upward is priority).
- `dwinmsgs` assumed to come from parent initialization process.
- assume `subfg` as a subgraph that can be modified by this function (add message factors)
  - should remove message prior factors from subgraph before returning.
- May modify `cliq` values.
  - `putMsgUpInit!(cliq, msg)`
  - `setCliqueStatus!(cliq, status)`
  - `setCliqDrawColor(cliq, "sienna")`
  - `notifyCliqDownInitStatus!(cliq, status)`

Algorithm:
- determine which downward messages influence initialization order
- initialize from singletons to most connected non-singletons
- revert back to needdownmsg if cycleInit does nothing
- can only ever return :initialized or :needdownmsg status

DevNotes
- TODO Lots of cleanup required, especially from calling function.
- TODO move directly into a CSM state function
"""
function doCliqInitDown!( subfg::AbstractDFG,
                          cliq::TreeClique,
                          initorder;
                          dbg::Bool=false,
                          logpath::String="/tmp/caesar/",
                          logger=ConsoleLogger() )
  #

  # store the cliqSubFg for later debugging
  if dbg
    DFG.saveDFG(subfg, joinpath(logpath,"logs/cliq$(cliq.index)/fg_beforedowninit"))
  end

  # cycle through vars and attempt init
  with_logger(logger) do
    @info "cliq $(cliq.index), doCliqInitDown! -- 5, cycle through vars and attempt init"
  end

  status = :needdownmsg
  if cycleInitByVarOrder!(subfg, initorder)
    status = :initialized
  end

  with_logger(logger) do
    @info "cliq $(cliq.index), doCliqInitDown! -- 6, current status: $status"
  end

  # store the cliqSubFg for later debugging
  if dbg
      DFG.saveDFG(subfg, joinpath(logpath,"logs/cliq$(cliq.index)/fg_afterdowninit"))
  end

  return status
end


function blockCliqUntilParentDownSolved(prnt::TreeClique; logger=ConsoleLogger())::Nothing
  @error("blockCliqUntilParentDownSolved is deprecated, use CSM directly")

  lbl = getLabel(prnt)

  with_logger(logger) do
    @info "blockCliqUntilParentDownSolved, prntcliq=$(prnt.index) | $lbl | going to fetch initdownchannel..."
  end
  flush(logger.stream)
  blockMsgDwnUntilStatus(prnt, :downsolved)

  return nothing
end

function blockMsgDwnUntilStatus(cliq::TreeClique, status::CliqStatus)
  @error("blockMsgDwnUntilStatus is deprecated, use CSM directly")
  while fetchMsgDwnInit(cliq).status != status
    wait(getSolveCondition(cliq))
  end
  nothing
end

# """
#    $SIGNATURES
#
# Determine if this `cliq` has been fully initialized and child cliques have completed their full upward inference.
# """
# function isCliqReadyInferenceUp(fgl::FactorGraph, tree::AbstractBayesTree, cliq::TreeClique)
#   isallinit = areCliqVariablesAllInitialized(fgl, cliq)
#
#   # check that all child cliques have also completed full up inference.
#   for chl in getChildren(tree, cliq)
#     isallinit &= isUpInferenceComplete(chl)
#   end
#   return isallinit
# end

# """
#     $SIGNATURES
#
# Perform cliq initalization calculation based on current state of the tree and factor graph,
# using upward message passing logic.
#
# Notes
# - adds msg priors added to clique subgraph
# - Return either of (:initialized, :upsolved, :needdownmsg, :badinit)
# - must use factors in cliq only, ensured by using subgraph -- TODO general case.
#
# DevNotes
# - FIXME, integrate with `8f. prepInitUp_StateMachine`
# """
# function doCliqAutoInitUpPart1!(subfg::AbstractDFG,
#                                 tree::AbstractBayesTree,
#                                 cliq::TreeClique;
#                                 up_solve_if_able::Bool=true,
#                                 multiproc::Bool=true,
#                                 logger=ConsoleLogger() )
#   #
#
#   # attempt initialize if necessary
#   if !areCliqVariablesAllInitialized(subfg, cliq)
#     # structure for all up message densities computed during this initialization procedure.
#     varorder = getCliqVarInitOrderUp(tree, cliq)
#     # do physical inits, ignore cycle return value
#     cycleInitByVarOrder!(subfg, varorder, logger=logger)
#   end
#
#   return nothing
# end


"""
    $SIGNATURES

Set the marginalized status of a clique.
"""
function setCliqAsMarginalized!(cliq::TreeClique, status::Bool)
  @warn "Busy deprecating setCliqAsMarginalized!, use CSM features directly"
  if status
    setCliqueStatus!(cliq, :marginalized)
  else
    if getCliqueData(cliq).initialized == :marginalized
      @info "Reverting clique $(cliq.index) to assumed :downsolved status"
      setCliqueStatus!(cliq, :downsolved)
    else
      error("Unknown clique de-marginalization requist for clique $(cliq.index), current status: $(cliq.initialized)")
    end
  end
end

"""
    $SIGNATURES

Run through entire tree and set cliques as marginalized if all clique variables are marginalized.

Notes:
- TODO can be made fully parallel, consider converting for use with `@threads` `for`.
"""
function updateTreeCliquesAsMarginalizedFromVars!(fgl::AbstractDFG, tree::AbstractBayesTree)::Nothing
  @warn "Busy deprecating updateTreeCliquesAsMarginalizedFromVars!, use CSM features directly"
  for (clid, cliq) in getCliques(tree)
    if isCliqMarginalizedFromVars(fgl, cliq)
      setCliqAsMarginalized!(cliq, true)
    end
  end
  nothing
end


##==============================================================================
## Delete at end v0.15.x
##==============================================================================


@deprecate putCliqueMsgDown!(cdata::BayesTreeNodeData, msg::LikelihoodMessage; from::Symbol=:nothing) putMsgDwnThis!(cdata, msg, from=from)

# """
#     $SIGNATURES

# Initialization downward message passing is different from regular inference since
# it is possible that none of the child cliq variables have been initialized.

# Notes
# - init upward msgs are individually stored in child cliques (pull model good).
# - fresh product of overlapping beliefs are calculated on each function call.
# - Assumed that `prnt` of siblings

# Dev Notes
# - This should be the initialization cycle of parent, build up bit by bit...
# """
# function prepCliqInitMsgsDown!(fgl::AbstractDFG,
#                                tree::AbstractBayesTree,
#                                prnt::TreeClique,
#                                cliq::TreeClique;
#                                logger=ConsoleLogger() )
#   #
#   # tt = split(string(now()), 'T')[end]
#   # with_logger(logger) do
#   #   @info "$(tt) prnt $(prnt.index), prepCliqInitMsgsDown! -- with cliq $(cliq.index)"
#   # end
  
#   # FIXME use only LikelihoodMessage
#   # check if any msgs should be multiplied together for the same variable
#   # msgspervar = LikelihoodMessage() # or maybe Dict{Int, LikelihoodMessage}()
#   msgspervar = getMsgInitDwnParent(tree, cliq, logger=logger)
#   # reference to default dict location
#   #JT 459 products = getMsgDwnThisInit(prnt)
#   products = getfetchCliqueInitMsgDown(prnt.data, from=:getMsgDwnThisInit) |> deepcopy
  
#   ## TODO use parent factors too
#   # intersect with the asking clique's separator variables
#   condenseDownMsgsProductPrntFactors!(fgl, products, msgspervar, prnt, cliq, logger)
  
#   # with_logger(logger) do
#   #   @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- vars fw/ down msgs=$(collect(keys(msgspervar)))"
#   # end
#   # flush(logger.stream)

#   # remove msgs that have no data
#   rmlist = Symbol[]
#   for (prsym,belmsg) in products.belief
#     if belmsg.inferdim < 1e-10
#       # no information so remove
#       push!(rmlist, prsym)
#     end
#   end
#   with_logger(logger) do
#     @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- rmlist, no inferdim, keys=$(rmlist)"
#   end
#   for pr in rmlist
#     delete!(products.belief, pr)
#   end

#   with_logger(logger) do
#     @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- product keys=$(collect(keys(products.belief)))"
#   end

#   # now put the newly computed message in the appropriate container
#   # FIXME THIS IS A PUSH MODEL, see #674 -- make pull model
#   putCliqueInitMsgDown!(getCliqueData(prnt), products)

#   return products
# end


# used during nonparametric CK preparation, when information from multiple siblings must be shared together
# const IntermediateSiblingMessagesTB{T} = Vector{TreeBelief{T}}
# const IntermediateMultiSiblingMessagesTB{T} = Dict{Symbol, Vector{TreeBelief{T}}}


# FIXME, better standardize intermediate types
# can be replaced by Vector{TreeBelief}
# const IntermediateSiblingMessages = Vector{Tuple{BallTreeDensity,Float64}}
# const IntermediateMultiSiblingMessages = Dict{Symbol, IntermediateSiblingMessages}


# Helper function for prepCliqInitMsgsDown!
# populate products with products of upward messages
function condenseDownMsgsProductOnly!(fgl::AbstractDFG,
                                      products::LikelihoodMessage,
                                      msgspervar::Dict{Symbol, <:AbstractVector}  )
  #
  error("condenseDownMsgsProductOnly!(::AbstractDFG,::LikelihoodMessage, ::IntermediateMultiSiblingMessages) is obsolete")
  # multiply multiple messages together
  for (msgsym, msgsBo) in msgspervar
    # check if this particular down message requires msgsym
    if exists(fgl, msgsym) # DFG.hasVariable(fgl, msgsym)
      if length(msgspervar[msgsym]) > 1
        msgs = getindex.(msgsBo, 1)
        haspars = 0.0
        for mb in msgsBo, val in mb[2]
          haspars += val
        end
        products[msgsym] = (manifoldProduct(msgs, getManifolds(fgl, msgsym)), haspars)
      else
        # transfer if only have a single belief
        products[msgsym] = (msgsBo[1][1], msgsBo[1][2])
      end
    else
      # not required, therefore remove from message to avoid confusion
      if haskey(products, msgsym)
        delete!(products, msgsym)
      end
    end
  end
  nothing
end

@deprecate getfetchCliqueMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing) getfetchCliqueInitMsgDown(cdata, from=from)

export getIdx

"""
    $SIGNATURES

Return interger index of desired variable element.

Example
-------
```julia
pp = RoME.Point2()
getIdx(pp, :posY) # = 2
```

Internal Notes
--------------
- uses number i < 100 for index number, and
- uses +100 offsets to track the minibatch number of the requested dimension
"""
function getIdx(pp::Tuple,
                sym::Symbol,
                i::Int=0)
  #
  error("getIdx is obsolete, use DistributedFactorGraphs objects/methods instead.")
  i-=100
  for p in pp
    i,j = getIdx(p, sym, i)
    if i > 0
      return i, j
    end
  end
  return i,-1
end
getIdx(pp::Symbol, sym::Symbol, i::Int=0) = pp==sym ? (abs(i)%100+1, div(abs(i)-100,100)) : (i-1, div(abs(i)-100,100))
function getIdx(pp::InferenceVariable, sym::Symbol, i::Int=0)
  return getIdx(pp.dimtype, sym)
end


# """
#     $SIGNATURES

# Get the latest down message from the parent node (without calculating anything).

# Notes
# - Different from down initialization messages that do calculate new values -- see `prepCliqInitMsgsDown!`.
# - Basically converts function `getDwnMsgs` from `Dict{Symbol,BallTreeDensity}` to `Dict{Symbol,Vector{BallTreeDensity}}`.

# Related

# getMsgInitDwnParent
# """
# function getMsgDwnParent(treel::AbstractBayesTree, cliq::TreeClique)
#   downmsgs = IntermediateMultiSiblingMessages()
#   prnts = getParent(treel, cliq)
#   if 0 < length(prnts)
#     prnt = prnts[1]
#     prntmsgs = getDwnMsgs(prnt)
#     for (key, bel) in prntmsgs
#       if !haskey(downmsgs, key)
#         downmsgs[key] = IntermediateSiblingMessages()
#       end
#       # TODO insert true inferred dim
#       push!(downmsgs[key], bel)
#     end
#   end

#   return downmsgs
# end


# function getMsgDwnThisInit(cdat::BayesTreeNodeData) 
#   iifdepwarn("#459 replace with getfetchCliqueMsgDown", :getMsgDwnThisInit)
#   return cdat.downInitMsg
# end

@deprecate fetchDataElement(dfg::AbstractDFG, varsym::Symbol, lbl::Symbol) fetchDataJSON(dfg, varsym, lbl)


function addMsgFactors!(subfg::AbstractDFG,
                        msgs::Dict{Symbol, <:AbstractVector} )  #Dict{Symbol, Vector{Tuple{BallTreeDensity, Float64}}} )
  # msgs::
  # add messages as priors to this sub factor graph
  @warn "Tuple{KDE,Floa64} specific version of addMsgFactors! is deprecated, use LikelihoodMessage version instead."
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, dms) in msgs
    for treebelief in dms
      if msym in svars
        # TODO should be on manifold prior, not just generic euclidean prior -- okay since variable on manifold, but not for long term
        fc = addFactor!(subfg, [msym], MsgPrior(manikde!(treebelief), treebelief.inferdim), graphinit=false)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end
