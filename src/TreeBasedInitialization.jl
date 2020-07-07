


function isCliqInitialized(cliq::TreeClique)::Bool
  return getCliqueData(cliq).initialized in [:initialized; :upsolved]
end

function isCliqUpSolved(cliq::TreeClique)::Bool
  return getCliqueData(cliq).initialized == :upsolved
end



"""
    $SIGNATURES

Return the most likely  ordering for initializing factor (assuming up solve
sequence).

Notes:
- sorts id for increasing number of connected factors.
"""
function getCliqVarInitOrderUp(tree::BayesTree, cliq::TreeClique)
  # rules to explore dimension from one to the other?

  # get all variable ids and number of associated factors
  allids = getCliqAllVarIds(cliq)
  nfcts = getCliqNumAssocFactorsPerVar(cliq)

  # get priors and singleton message variables (without partials)
  prids = getCliqVarIdsPriors(cliq, getCliqAllVarIds(cliq), false)

  # get current up msgs in the init process (now have all singletons)
  upmsgs = getMsgsUpInitChildren(tree, cliq, TreeBelief)
  upmsgids = collect(keys(upmsgs))

  # all singleton variables
  singids = union(prids, upmsgids)

  # add msg marginal prior (singletons) to number of factors
  for msid in upmsgids
    nfcts[msid .== allids] .+= 1
  end

  # sort permutation order for increasing number of factor association
  nfctsp = sortperm(nfcts)
  sortedids = allids[nfctsp]

  # organize the prior variables separately with asceding factor count
  initorder = Symbol[]
  for id in sortedids
    if id in singids
      push!(initorder, id)
    end
  end
  # in ascending order of number of factors
  for id in sortedids
    if !(id in initorder)
      push!(initorder, id)
    end
  end
  return initorder
end


"""
    $SIGNATURES

Return true if clique has completed the local upward direction inference procedure.
"""
isUpInferenceComplete(cliq::TreeClique) = getCliqueData(cliq).upsolved

function areCliqVariablesAllInitialized(dfg::G, cliq::TreeClique) where {G <: AbstractDFG}
  allids = getCliqAllVarIds(cliq)
  isallinit = true
  for vid in allids
    var = DFG.getVariable(dfg, vid)
    isallinit &= isInitialized(var)
  end
  isallinit
end



"""
    $SIGNATURES

Return true if all variables in clique are considered marginalized (and initialized).
"""
function areCliqVariablesAllMarginalized(subfg::AbstractDFG,
                                         cliq::TreeClique)
  for vsym in getCliqAllVarIds(cliq)
    vert = getVariable(subfg, vsym)
    if !isMarginalized(vert) || !isInitialized(vert)
      return false
    end
  end
  return true
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
      # need to set the upward messages
      msgs = prepCliqInitMsgsUp(dfg, cliq)
      putMsgUpThis!(cliq, msgs)

      prnt = getParent(tree, cliq)
      if length(prnt) > 0
        # THIS IS FOR INIT PASSES ONLY
        putMsgUpInit!(cliq, msgs, logger)
      end

      setCliqStatus!(cliq, :marginalized)
      setCliqDrawColor(cliq, "blue")

      # set flag, looks to be previously unused???
      getCliqueData(cliq).allmarginalized = true
    end
  end
  nothing
end




"""
    $SIGNATURES

Cycle through var order and initialize variables as possible in `subfg::AbstractDFG`.
Return true if something was updated.

Notes:
- assumed `subfg` is a subgraph containing only the factors that can be used.
  - including the required up or down messages
- intended for both up and down initialization operations.

Dev Notes
- Should monitor updates based on the number of inferred & solvable dimensions
"""
function cycleInitByVarOrder!(subfg::AbstractDFG,
                              varorder::Vector{Symbol};
                              logger=ConsoleLogger()  )::Bool
  #
  with_logger(logger) do
    @info "cycleInitByVarOrder! -- varorder=$(varorder)"
  end
  retval = false
  count = 1
  while count > 0
    count = 0
    for vsym in varorder
      var = DFG.getVariable(subfg, vsym)
      isinit = isInitialized(var)
      with_logger(logger) do
        @info "var.label=$(var.label) is initialized=$(isinit)"
      end
      doautoinit!(subfg, [var;], logger=logger)
      if isinit != isInitialized(var)
        count += 1
        retval = true
      end
    end
  end
  with_logger(logger) do
    @info "cycleInitByVarOrder!, retval=$(retval)"
  end
  flush(logger.stream)
  return retval
end


# currently for internal use only
# initialize variables based on best current achievable ordering
# OBVIOUSLY a lot of refactoring and consolidation needed with cliqGibbs / approxCliqMarginalUp
function initSolveSubFg!(subfg::AbstractDFG,
                         logger=ConsoleLogger() )
  #
  varorder = getSubFgPriorityInitOrder(subfg, logger)
  with_logger(logger) do
    @info "initSolveSubFg! -- varorder=$varorder"
  end
  cycleInitByVarOrder!(subfg, varorder, logger=logger)
  nothing
end


# Helper function for prepCliqInitMsgsDown!
# populate products with products of upward messages
function condenseDownMsgsProductOnly!(fgl::AbstractDFG,
                                      products::LikelihoodMessage,
                                      msgspervar::IntermediateMultiSiblingMessages  )
  #
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



# Helper function for prepCliqInitMsgsDown!
# future, be used in a cached system with parent in one location only for all siblings
function condenseDownMsgsProductPrntFactors!(fgl::G,
                                             products::LikelihoodMessage,
                                             msgspervar::IntermediateMultiSiblingMessages,
                                             prnt::TreeClique,
                                             cliq::TreeClique,
                                             logger=ConsoleLogger()  ) where G <: AbstractDFG
  #

  # determine the variables of interest
  reqMsgIds = collect(keys(msgspervar))
  # unique frontals per cliq
  prntvars = intersect(getCliqSeparatorVarIds(cliq), getCliqAllVarIds(prnt))
  lvarids = union(prntvars, reqMsgIds)
  # determine allowable factors, if any (only from parent cliq)
  awfcts = getCliqFactorIdsAll(prnt)
  with_logger(logger) do
      @info "condenseDownMsgsProductPrntFactors! -- reqMsgIds=$(reqMsgIds),"
      @info "condenseDownMsgsProductPrntFactors! -- vars=$(lvarids),"
      @info "condenseDownMsgsProductPrntFactors! -- allow factors $awfcts"
  end

  # build required subgraph for parent/sibling down msgs
  lsfg = buildSubgraph(fgl, lvarids, 1)

  tempfcts = lsf(lsfg)
  dellist = setdiff(awfcts, tempfcts)
  for delf in dellist
    # TODO -- double check this deletefactor method is leaving the right parent sharing factor graph behind
    if exists(lsfg, delf) # hasFactor
      deleteFactor!(lsfg,delf)
    end
  end
  with_logger(logger) do
      @info "condenseDownMsgsProductPrntFactors! -- lsfg fcts=$(lsf(lsfg)),"
      @info "condenseDownMsgsProductPrntFactors! -- excess factors $dellist"
  end

  # add message priors
  addMsgFactors!(lsfg, msgspervar)

  # perform initialization/inference
  # ....uhhh TODO
  initSolveSubFg!(lsfg, logger)

      # QUICK DBG CODE
      vars = ls(lsfg)
      len = length(vars)
      tdims = Vector{Float64}(undef, len)
      isinit = Vector{Bool}(undef, len)
      for idx in 1:len
        tdims[idx] = getVariableSolvableDim(lsfg, vars[idx])
        isinit[idx] = isInitialized(lsfg, vars[idx])
      end
      with_logger(logger) do
          @info "condenseDownMsgsProductPrntFactors! -- after cycle init: vars=$vars"
          @info "condenseDownMsgsProductPrntFactors! -- after cycle init: tdims=$tdims"
          @info "condenseDownMsgsProductPrntFactors! -- after cycle init: isinit=$isinit"
      end
      # QUICK DBG CODE

  # extract complete downward marginal msg priors
  for id in intersect(getCliqSeparatorVarIds(cliq), lvarids)
    vari = getVariable(lsfg, id)
    products.belief[id] = TreeBelief(vari)
    # products.belief[id] = (getKDE(vari), getVariableInferredDim(vari))
  end

  # if recycling lsfg
  # deleteMsgFactors!(...)

  nothing
end

"""
    $SIGNATURES

Initialization downward message passing is different from regular inference since
it is possible that none of the child cliq variables have been initialized.

Notes
- init upward msgs are individually stored in child cliques ().
- fresh product of overlapping beliefs are calculated on each function call.
- Assumed that `prnt` of siblings

Dev Notes
- This should be the initialization cycle of parent, build up bit by bit...
"""
function prepCliqInitMsgsDown!(fgl::AbstractDFG,
                               tree::AbstractBayesTree,
                               prnt::TreeClique,
                               cliq::TreeClique;
                               logger=ConsoleLogger(),
                               dbgnew::Bool=true  )
  #
  tt = split(string(now()), 'T')[end]
  with_logger(logger) do
    @info "$(tt) prnt $(prnt.index), prepCliqInitMsgsDown! -- with cliq $(cliq.index)"
  end
  # get the current messages ~~stored in~~ [going to] the parent
  currmsgs = getMsgsUpInitChildren(tree, prnt, TreeBelief, [cliq.index;])
  with_logger(logger) do
    @info "prnt $(prnt.index), prepCliqInitMsgsDown! -- msg ids::Int=$(collect(keys(currmsgs)))"
  end

  # FIXME drop IntermediateMultiSiblingMessages and use only LikelihoodMessage
  # check if any msgs should be multiplied together for the same variable
  # msgspervar = LikelihoodMessage() # or maybe Dict{Int, LikelihoodMessage}()
  msgspervar = IntermediateMultiSiblingMessages()

  for (msgcliqid, msgs) in currmsgs
    with_logger(logger) do
      @info "prepCliqInitMsgsDown! -- msgcliqid=$msgcliqid, msgs.belief=$(collect(keys(msgs.belief)))"
    end
    for (msgsym, msg) in msgs.belief
      if !haskey(msgspervar, msgsym)
        # there will be an entire list...
        msgspervar[msgsym] = IntermediateSiblingMessages()
      end
      with_logger(logger) do
        @info "prepCliqInitMsgsDown! -- msgcliqid=$(msgcliqid), msgsym $(msgsym), inferdim=$(msg.inferdim)"
      end
      push!(msgspervar[msgsym], msg)
    end
  end

  with_logger(logger) do
    @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- vars fw/ down msgs=$(collect(keys(msgspervar)))"
  end

  flush(logger.stream)

  # reference to default dict location
  products = getMsgDwnThisInit(prnt)

  ## TODO use parent factors too
  # intersect with the asking clique's separator variables

    # products only method
    # @assert dbgnew "must use dbgnew=true for tree down init msg"
    if dbgnew
        condenseDownMsgsProductPrntFactors!(fgl, products, msgspervar, prnt, cliq, logger) # WIP -- not ready yet
    else
        condenseDownMsgsProductOnly!(fgl, products, msgspervar) # BASELINE deprecated
    end

  # remove msgs that have no data
  rmlist = Symbol[]
  for (prsym,belmsg) in products.belief
    if belmsg.inferdim < 1e-10
      # no information so remove
      push!(rmlist, prsym)
    end
  end
  with_logger(logger) do
    @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- rmlist, no inferdim, keys=$(rmlist)"
  end
  for pr in rmlist
    delete!(products.belief, pr)
  end

  with_logger(logger) do
    @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- product keys=$(collect(keys(products.belief)))"
  end
  return products
end


"""
    $SIGNATURES

Prepare the upward inference messages from clique to parent and return as `Dict{Symbol}`.

Notes
- Does not require tree message likelihood factors in subfg.
- Also see #579 regarding elimited likelihoods and priors.
"""
function prepCliqInitMsgsUp(subfg::AbstractDFG,
                            cliq::TreeClique,
                            logger=ConsoleLogger() )::LikelihoodMessage
  #
  # construct init's up msg to place in parent from initialized separator variables
  msg = LikelihoodMessage()
  seps = getCliqSeparatorVarIds(cliq)
  with_logger(logger) do
    @info "prepCliqInitMsgsUp, seps=$seps"
  end
  for vid in seps
    var = DFG.getVariable(subfg, vid)
    if isInitialized(var)
      msg.belief[Symbol(var.label)] = TreeBelief(var)
    end
  end
  return msg
end



function blockCliqUntilParentDownSolved(prnt::TreeClique; logger=ConsoleLogger())::Nothing
  #
  lbl = getLabel(prnt)

  with_logger(logger) do
    @info "blockCliqUntilParentDownSolved, prntcliq=$(prnt.index) | $lbl | going to fetch initdownchannel..."
  end
  flush(logger.stream)
  blockMsgDwnUntilStatus(prnt, :downsolved)

  return nothing
end



"""
    $SIGNATURES

Block the thread until child cliques of `prnt::TreeClique` have finished
attempting upward initialization -- i.e. have status result.
Return `::Dict{Symbol}` indicating whether next action that should be taken
for each child clique.

Notes:
- See status options at `getCliqStatusUp(..)`.
- Can be called multiple times
"""
function blockCliqUntilChildrenHaveUpStatus(tree::AbstractBayesTree,
                                            prnt::TreeClique,
                                            logger=ConsoleLogger() )
  #
  ret = Dict{Int, Symbol}()
  chlr = getChildren(tree, prnt)
  for ch in chlr
    # either wait to fetch new result, or report or result
    chst = getCliqStatusUp(ch)
    with_logger(logger) do
      @info "cliq $(prnt.index), child $(ch.index) status is $(chst), isready(initUpCh)=$(isready(getMsgUpInitChannel_(ch)))."
    end
    flush(logger.stream)
    ret[ch.index] = fetchMsgUpInit(ch).status
  end
  with_logger(logger) do
      tt = split(string(now()), 'T')[end]
      @info "cliq $(prnt.index), $(tt), fetched all, keys=$(keys(ret))."
  end
  return ret
end

"""
    $SIGNATURES

Wait here if all siblings and the parent status are `:needdownmsg`.
Return true when parent is `:initialized` after all were `:needdownmsg`

Notes
- used for regulating long need down message chains.
- exit strategy is parent becomes status `:initialized`.
"""
function blockCliqSiblingsParentNeedDown(tree::AbstractBayesTree,
                                         cliq::TreeClique; logger=ConsoleLogger())
  #
  with_logger(logger) do
    @info "cliq $(cliq.index), blockCliqSiblingsParentNeedDown -- start of function"
  end
  # ret = Dict{Int, Symbol}()
  # flush(logger.stream)
  prnt = getParent(tree, cliq)
  allneeddwn = true
  if length(prnt) > 0
    prstat = getCliqStatus(prnt[1])
    if prstat == :needdownmsg
      for ch in getChildren(tree, prnt[1])
        chst = getCliqStatusUp(ch)
        if chst != :needdownmsg
          allneeddwn = false
          break;
        end
      end

      if allneeddwn
        with_logger(logger) do
          tt = split(string(now()), 'T')[end]
          @warn "$(tt) | $(current_task()), cliq $(cliq.index), block since all siblings/parent($(prnt[1].index)) :needdownmsg."
        end
        flush(logger.stream)
        # do actual fetch
        prtmsg = fetchMsgDwnInit(prnt[1]).status
        with_logger(logger) do
            tt = split(string(now()), 'T')[end]
          @info "$tt | $(current_task()) clique $(prnt[1].index), blockCliqSiblingsParentNeedDown -- after fetch $prstat, $prtmsg"
        end
        if prtmsg == :initialized
          return true
        else
            with_logger(logger) do
                tt = split(string(now()), 'T')[end]
                @warn "$tt | $(current_task()) Clique $(prnt[1].index), maybe clear down init message $prtmsg"
            end
        end
      end
    end
  end
  return false
end



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
# - FIXME, integrate with `8f. mustInitUpCliq_StateMachine`
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

function printCliqInitPartialInfo(subfg, cliq, logger=ConsoleLogger())
  varids = getCliqAllVarIds(cliq)
  initstatus = Vector{Bool}(undef, length(varids))
  initpartial = Vector{Float64}(undef, length(varids))
  for i in 1:length(varids)
    initstatus[i] = getSolverData(getVariable(subfg, varids[i])).initialized
    initpartial[i] = getSolverData(getVariable(subfg, varids[i])).inferdim
  end
  with_logger(logger) do
    tt = split(string(now()),'T')[end]
    @info "$tt, cliq $(cliq.index), PARINIT: $varids | $initstatus | $initpartial"
  end
end

"""
    $SIGNATURES

Notify of new up status and message.

Notes
- Major part of #459 consolidation effort.
"""
function prepPutCliqueStatusMsgUp!(csmc::CliqStateMachineContainer,
                                   status::Symbol  )
  #
  # construct init's up msg from initialized separator variables
  upinitmsg = prepCliqInitMsgsUp(csmc.cliqSubFg, csmc.cliq)
  # put the init upinitmsg
  putMsgUpInit!(csmc.cliq, upinitmsg, csmc.logger)
  if getCliqStatus(csmc.cliq) != status
	infocsm(csmc, "prepPutCliqueStatusMsgUp! -- notify status=$status")
	notifyCliqUpInitStatus!(csmc.cliq, status, logger=csmc.logger)
  end

  # print a little late
  with_logger(csmc.logger) do
    tt = split(string(now()),'T')[end]
    @info "$tt, cliq $(csmc.cliq.index), 8g, doCliqUpsSolveInit. -- postupinitmsg with $(collect(keys(upinitmsg.belief)))"
  end

  # return new up messages in case the user wants to see
  return upinitmsg
end

"""
    $SIGNATURES

Follows cliq initalization calculation and attempts full upsolve
based on current state of the tree and factor graph,
using upward message passing logic.

Notes
- Removes msg priors added to clique subgraph
- Return either of (:initialized, :upsolved, :needdownmsg, :badinit)
- must use factors in cliq only, ensured by using subgraph -- TODO general case.
"""
function doCliqAutoInitUpPart2!(csmc::CliqStateMachineContainer;
                                up_solve_if_able::Bool=true,
                                multiproc::Bool=true,
                                logger=ConsoleLogger()  )
  #

  cliqst = getCliqStatus(csmc.cliq)
  status = (cliqst == :initialized || length(getParent(csmc.tree, csmc.cliq)) == 0) ? cliqst : :needdownmsg

  # check if all cliq vars have been initialized so that full inference can occur on clique
  if areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq)
    with_logger(logger) do
      tt = split(string(now()),'T')[end]
      @info "$(tt), cliq $(csmc.cliq.index), doCliqUpSolvePart2!, clique status = $(status)"
    end

    # TODO replace with msg channels only
    urt = approxCliqMarginalUp!(csmc, logger=csmc.logger)
    # is clique fully upsolved or only partially?
    # TODO verify the need for this update (likely part of larger refactor, WIP #459)
    putMsgUpThis!(csmc.cliq, urt.keepupmsgs)
    updateFGBT!(csmc.cliqSubFg, csmc.cliq, urt, dbg=getSolverParams(csmc.cliqSubFg).dbg, fillcolor="brown", logger=csmc.logger)

    # set clique color accordingly, using local memory
    setCliqDrawColor(csmc.cliq, isCliqFullDim(csmc.cliqSubFg, csmc.cliq) ? "pink" : "tomato1")
    getCliqueData(csmc.cliq).upsolved = true
    status = :upsolved

    # NOTE moving message update here before moving back to CSM.jl

  else
    with_logger(logger) do
      @info "cliq $(csmc.cliq.index), all variables not initialized, status = $(status)"
    end
  end

  return status
end


"""
    $SIGNATURES

Special function to do initialization in downward direction, assuming that not all
variables can be initialized.  Relies on outside down messages.

Notes:
- assumed this `cliq` is being initialized from a previous `:needdownmsg` status.
- will use all possible local factors of cliquq in initilization process
- similar to upward initialization, but uses different message structure
  - first draft assumes upward messages will not be used,
  - full up solve still required which explicitly depends on upward messages.

Dev Notes
- Streamline add/delete msg priors from calling function and csm.
- TODO replace with nested 'minimum degree' type variable ordering.
"""
function getCliqInitVarOrderDown(dfg::G,
                                 cliq::TreeClique,
                                 downmsgs::LikelihoodMessage )::Vector{Symbol} where G <: AbstractDFG
  #
  allsyms = getCliqAllVarIds(cliq)
  # convert input downmsg var symbols to integers (also assumed as prior beliefs)
  # make sure ids are in the clique set, since parent may have more variables.
  dwnmsgsym = intersect(collect(keys(downmsgs.belief)), DFG.listVariables(dfg)) #dfg.IDs
  # dwnmsgids =  map(x -> dfg.IDs[x], dwnmsgsym )
  dwnvarids = intersect(allsyms, dwnmsgsym)

  # find any other prior factors (might have partials)
  prvarids = getCliqVarIdsPriors(cliq, allsyms, true)
  hassinglids = union(dwnvarids, prvarids)

  # Get all other variable factor counts
  nfcts = getCliqNumAssocFactorsPerVar(cliq)
  # add msg marginal prior (singletons) to number of factors
  for msid in dwnmsgsym
    nfcts[msid .== allsyms] .+= 1
  end

  # sort permutation order for increasing number of factor association
  nfctsp = sortperm(nfcts)
  sortedids = allsyms[nfctsp]

  # all singleton variables
  singids = union(prvarids, dwnvarids)

  # organize the prior variables separately with asceding factor count
  initorder = Symbol[] #zeros(Int, 0)
  for id in sortedids
    if id in singids
      push!(initorder, id)
    end
  end
  # sort remaining variables for increasing associated factors
  for id in sortedids
    if !(id in initorder)
      push!(initorder, id)
    end
  end

  # return variable order
  return initorder
end


"""
    $SIGNATURES

Return true or false depending on whether child cliques are all up solved.
"""
function areCliqChildrenAllUpSolved(treel::AbstractBayesTree,
                                    prnt::TreeClique)::Bool
  #
  for ch in getChildren(treel, prnt)
    if !isCliqUpSolved(ch)
      return false
    end
  end
  return true
end


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
  - `setCliqStatus!(cliq, status)`
  - `setCliqDrawColor(cliq, "sienna")`
  - `notifyCliqDownInitStatus!(cliq, status)`

Algorithm:
- determine which downward messages influence initialization order
- initialize from singletons to most connected non-singletons
- revert back to needdownmsg if cycleInit does nothing
- can only ever return :initialized or :needdownmsg status

DevNotes
- TODO Lots of cleanup required, especially from calling function.
"""
function doCliqInitDown!(subfg::AbstractDFG,
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


"""
    $SIGNATURES

Return `true` if any of the children cliques have status `:needdownmsg`.
"""
function areCliqChildrenNeedDownMsg(children::Vector{TreeClique})::Bool
  for ch in children
    if getCliqStatus(ch) == :needdownmsg
      return true
    end
  end
  return false
end

function areCliqChildrenNeedDownMsg(tree::AbstractBayesTree, cliq::TreeClique)::Bool
  areCliqChildrenNeedDownMsg( getChildren(tree, cliq) )
end


"""
    $SIGNATURES

Return true if has parent with status `:needdownmsg`.
"""
function isCliqParentNeedDownMsg(tree::AbstractBayesTree, cliq::TreeClique, logger=ConsoleLogger())
  prnt = getParent(tree, cliq)
  if length(prnt) == 0
    return false
  end
  prstat = getCliqStatus(prnt[1])
  with_logger(logger) do
    @info "$(current_task()) Clique $(cliq.index), isCliqParentNeedDownMsg -- parent status: $(prstat)"
  end
  return prstat == :needdownmsg
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
