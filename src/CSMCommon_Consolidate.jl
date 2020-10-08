
## ===================================================================================================================
##  CSM logging functions
## ===================================================================================================================

"""
    $SIGNATURES

Specialized info logger print function to show clique state machine information
in a standardized form.
"""
function infocsm(csmc::CliqStateMachineContainer, str::A) where {A <: AbstractString}

  tm = string(Dates.now())
  tmt = split(tm, 'T')[end]

  lbl = getLabel(csmc.cliq)
  lbl1 = split(lbl,',')[1]
  cliqst = getCliqueStatus(csmc.cliq)

  with_logger(csmc.logger) do
    @info "$tmt | $(csmc.cliq.index)---$lbl1 @ $(cliqst) | "*str
  end
  flush(csmc.logger.stream)
  nothing
end

"""
    $SIGNATURES
Helper function to log a message at a specific level to a clique identified by `csm_i` where i = cliq.index
Notes:
- Related to infocsm.
- Different approach to logging that uses the build in logging functionality to provide more flexibility.
- Can be used with LoggingExtras.jl 
"""
function logCSM(csmc, msg::String; loglevel::Logging.LogLevel=Logging.Debug, maxlog=nothing, kwargs...)
  #Debug = -1000
  #Info = 0
  #Warn = 1000
  #Error = 2000
  @logmsg(loglevel,
          msg,
          _module=begin
            bt = backtrace()
            funcsym=(:logCSM, Symbol("logCSM##kw")) #always use the calling function of logCSM
            frame, caller = Base.firstcaller(bt, funcsym)
            # TODO: Is it reasonable to attribute callers without linfo to Core?
            caller.linfo isa Core.MethodInstance ? caller.linfo.def.module : Core
          end,
          _file=String(caller.file),
          _line=caller.line,
          _id=(frame,funcsym),
          # caller=caller, 
          # st4 = stacktrace()[4],
          _group = Symbol("csm_$(csmc.cliq.index)"),
          maxlog=maxlog,
          kwargs...)
  
  return nothing
end

## ===================================================================================================================
## CSM Monitor functions
## ===================================================================================================================

"""
    $SIGNATURES
Monitor CSM tasks for failures and propagate error to the other CSMs to cleanly exit. 
"""
function monitorCSMs(tree, alltasks; forceIntExc::Bool=false)
  task = @async begin
    while true
      all(istaskdone.(alltasks)) && (@info "monitorCSMs: all tasks done"; break)
      for (i,t) in enumerate(alltasks)
        if istaskfailed(t)
          if forceIntExc
            @error "Task $i failed, sending InterruptExceptions to all running CSM tasks"
            throwIntExcToAllTasks(alltasks)
            @debug "done with throwIntExcToAllTasks"
          else
            @error "Task $i failed, sending error to all cliques"
            bruteForcePushErrorCSM(tree)
            # for tree.messages
            @info "All cliques should have exited"
          end
        end
      end
      sleep(1)
    end
  end
  return task
end

function throwIntExcToAllTasks(alltasks)
  for (i,t) in enumerate(alltasks)
    if !istaskdone(alltasks[i])
      @debug "Sending InterruptExceptions to CSM task $i"
      schedule(alltasks[i], InterruptException(), error=true)
      @debug "InterruptExceptions CSM task $i"
    end
  end 
  return nothing
end

function bruteForcePushErrorCSM(tree)
    errMsg = LikelihoodMessage(status=:ERROR_STATUS)
    for (i, ch) in tree.messages

        if isready(ch.upMsg)
            take!(ch.upMsg)
        else
            @info("Up edge $i", ch.upMsg)
            @async put!(ch.upMsg, errMsg)
        end
        if isready(ch.downMsg)
            take!(ch.downMsg)
        else
            @info("Down edge $i", ch.downMsg)
            @async put!(ch.downMsg, errMsg)
        end

    end

    for (i, ch) in tree.messages

        while isready(ch.upMsg)
            @info "cleanup take on $i up"
            take!(ch.upMsg)
        end
        while isready(ch.downMsg)
            @info "cleanup take on $i down"
            take!(ch.downMsg)
        end

    end

end


## ==================================================================================================================='
# Common nonparametric and parametric CSM functions
## ==================================================================================================================='



"""
    $SIGNATURES

One of the last steps in CSM to clean up after a down solve.

Notes
- CSM function 11b.
"""
function cleanupAfterDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  # RECENT split from 11 (using #760 solution for deleteMsgFactors)
  opts = getSolverParams(csmc.cliqSubFg)

  # set PPE and solved for all frontals
  for sym in getCliqFrontalVarIds(csmc.cliq)
    # set PPE in cliqSubFg
    setVariablePosteriorEstimates!(csmc.cliqSubFg, sym)
    # set solved flag
    vari = getVariable(csmc.cliqSubFg, sym)
    setSolvedCount!(vari, getSolvedCount(vari, :default)+1, :default )
  end

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterdownsolve")

  # transfer results to main factor graph
  frsyms = getCliqFrontalVarIds(csmc.cliq)
  infocsm(csmc, "11, finishingCliq -- going for transferUpdateSubGraph! on $frsyms")
  transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger, updatePPE=true)

  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- before notifyCliqDownInitStatus!")
  notifyCliqDownInitStatus!(csmc.cliq, :downsolved, logger=csmc.logger)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- just notified notifyCliqDownInitStatus!")

  # remove msg factors that were added to the subfg
  rmFcts = deleteMsgFactors!(csmc.cliqSubFg)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- removing all up/dwn message factors, length=$(length(rmFcts))")

  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- finished, exiting CSM on clique=$(csmc.cliq.index)")
  # and finished
  return IncrementalInference.exitStateMachine
end



"""
$SIGNATURES

Root clique upsolve and downsolve are equivalent, so skip a repeat downsolve, just set messages and just exit directly.

Notes
- State machine function nr. 10b
- Separate out during #459 dwnMsg consolidation.

DevNotes
- TODO should this consolidate some work with 11b?
"""
function specialCaseRootDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  # this is the root clique, so assume already downsolved -- only special case
  dwnmsgs = getCliqDownMsgsAfterDownSolve(csmc.cliqSubFg, csmc.cliq)
  setCliqDrawColor(csmc.cliq, "lightblue")

  # this part looks like a pull model
  # JT 459 putMsgDwnThis!(csmc.cliq, dwnmsgs)
  putDwnMsgConsolidated!(csmc.cliq.data, dwnmsgs) # , from=:putMsgDwnThis!  putCliqueMsgDown!
  setCliqueStatus!(csmc.cliq, :downsolved)
  csmc.dodownsolve = false

  # Update estimates and transfer back to the graph
  frsyms = getCliqFrontalVarIds(csmc.cliq)
  # set PPE and solved for all frontals
  for sym in frsyms
    # set PPE in cliqSubFg
    setVariablePosteriorEstimates!(csmc.cliqSubFg, sym)
    # set solved flag
    vari = getVariable(csmc.cliqSubFg, sym)
    setSolvedCount!(vari, getSolvedCount(vari, :default)+1, :default )
  end
  # Transfer to parent graph
  transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, updatePPE=true)

  prepPutCliqueStatusMsgDwn!(csmc, :downsolved)
  # notifyCliqDownInitStatus!(csmc.cliq, :downsolved, logger=csmc.logger)

  # bye
  return IncrementalInference.exitStateMachine
end


"""
$SIGNATURES

Quick redirection of out-marginalized cliques to downsolve path, or wait on children cliques to get a csm status.

Notes
- State machine function nr.4
"""
function canCliqMargSkipUpSolve_StateMachine(csmc::CliqStateMachineContainer)

  cliqst = getCliqueStatus(csmc.oldcliqdata)
  infocsm(csmc, "4, canCliqMargSkipUpSolve_StateMachine, $cliqst, csmc.incremental=$(csmc.incremental)")

  # if clique is out-marginalized, then no reason to continue with upsolve
  # marginalized state is set in `canCliqMargRecycle_StateMachine`
  if cliqst == :marginalized
    # go to 10 -- Add case for IIF issue #474
    return canCliqDownSolve_StateMachine
  end

  # go to 4e
  return blockUntilChildrenHaveStatus_StateMachine
end



"""
    $SIGNATURES

Placeholder function part of #459 dwnMsg consolidation to send current up message, part of :needdownmsg downinit cascading.

Notes
- State machine function nr. 8k
"""
function sendCurrentUpMsg_StateMachine(csmc::CliqStateMachineContainer)
  # set messages if children :needdownmsg
  infocsm(csmc, "8k, sendCurrentUpMsg_StateMachine -- must set messages for future down init")
  # construct init's up msg to place in parent from initialized separator variables

  # consolidated up messaging (#459)
  infocsm(csmc, "8k, sendCurrentUpMsg_StateMachine -- putting fake upinitmsg in this cliq")
  upmsg = prepCliqueMsgUpConsolidated(csmc.cliqSubFg, csmc.cliq, getCliqueStatus(csmc.cliq), logger=csmc.logger)
  prepPutCliqueStatusMsgUp!(csmc, upmsg=upmsg)

  # also send a down message -- seem weird while doing #459 but okay
  cliqst = getCliqueStatus(csmc.cliq)
  notifyCliqDownInitStatus!(csmc.cliq, cliqst, logger=csmc.logger)

  # Legend: initialized but not solved yet (likely child cliques that depend on downward autoinit msgs),
  setCliqDrawColor(csmc.cliq, "sienna")

  infocsm(csmc, "8k, sendCurrentUpMsg_StateMachine -- near-end down init attempt, $cliqst.")

  # go to 8b
  return attemptCliqInitUp_StateMachine
end



"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2
"""
function buildCliqSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  infocsm(csmc, "2, build subgraph syms=$(getCliqAllVarIds(csmc.cliq))")
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, csmc.cliq)

  # if dfg, store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_build")

  # go to 4
  return canCliqMargSkipUpSolve_StateMachine
end


"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2r
"""
function buildCliqSubgraphForDown_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarIds(csmc.cliq)
  infocsm(csmc, "2r, build subgraph syms=$(syms)")
  csmc.cliqSubFg = buildSubgraph(csmc.dfg, syms, 1; verbose=false)

  opts = getSolverParams(csmc.dfg)
  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_build_down")

  # go to 10
  return canCliqDownSolve_StateMachine
end

"""
    $SIGNATURES

Either construct and notify of a new upward initialization message and progress to downsolve checks,
or circle back and start building the local clique subgraph.

Notes
- State machine function nr.1
- Root clique message should be empty since it has an empty separator.
"""
function isCliqUpSolved_StateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "1, isCliqUpSolved_StateMachine")
  cliqst = getCliqueStatus(csmc.cliq)

  # if upward complete for any reason, prepare and send new upward message
  if cliqst in [:upsolved; :downsolved; :marginalized; :uprecycled]
    # construct init's up msg from initialized separator variables
    # NOTE cliqSubFg has not been copied yet
    prepPutCliqueStatusMsgUp!(csmc, cliqst, dfg=csmc.dfg)
    #go to 10
    return canCliqDownSolve_StateMachine
  end
  # go to 2
  return buildCliqSubgraph_StateMachine
end


"""
    $SIGNATURES

Final determination on whether can promote clique to `:uprecycled`.

Notes
- State machine function nr.0b
- Assume children clique status is available
- Will return to regular init-solve if new information in children -- ie not uprecycle or marginalized
"""
function checkChildrenAllUpRecycled_StateMachine(csmc::CliqStateMachineContainer)
  count = Int[]
  chldr = getChildren(csmc.tree, csmc.cliq)
  for ch in chldr
    chst = getCliqueStatus(ch)
    if chst in [:uprecycled; :marginalized]
      push!(count, 1)
    end
  end
  infocsm(csmc, "0b, checkChildrenAllUpRecycled_StateMachine -- length(chldr)=$(length(chldr)), sum(count)=$(sum(count))")

  # all children can be used for uprecycled -- i.e. no children have new information
  if sum(count) == length(chldr)
    # set up msg and exit go to 1
    sdims = Dict{Symbol,Float64}()
    for varid in getCliqAllVarIds(csmc.cliq)
      sdims[varid] = 0.0
    end

    # NOTE busy consolidating #459
    updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
        # setCliqueStatus!(csmc.cliq, :uprecycled)
        # replacing similar functionality from CSM 1.
    if getSolverParams(csmc.cliqSubFg).dbg
      tmnow = now()
      tmpst = getCliqueStatus(csmc.cliq)
      @async begin
        mkpath(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)"))
        open(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)","incremental.log"), "w") do f
          println(f, tmnow, ", marginalized from previous status ", tmpst)
        end
      end
    end
    prepPutCliqueStatusMsgUp!(csmc, :uprecycled, dfg=csmc.dfg)
    setCliqDrawColor(csmc.cliq, "orange")
    #go to 10
    return canCliqDownSolve_StateMachine
        # # go to 1
        # return isCliqUpSolved_StateMachine
  end

  # return to regular solve, go to 2
  return buildCliqSubgraph_StateMachine
end

"""
    $SIGNATURES

Determine if clique is upsolved by incremental update and exit the state machine.

Notes
- State machine function nr.0c
- can recycle if two checks:
  - previous clique was identically downsolved
  - all children are also :uprecycled
"""
function canCliqIncrRecycle_StateMachine(csmc::CliqStateMachineContainer)
	# check if should be trying and can recycle clique computations
    if csmc.incremental && getCliqueStatus(csmc.oldcliqdata) == :downsolved
      csmc.cliq.data.isCliqReused = true
      # check if a subgraph will be needed later
      if csmc.dodownsolve
        # yes need subgraph and need more checks, so go to 2
        return buildCliqSubgraph_StateMachine
      else
        # one or two checks say yes, so go to 4
        return canCliqMargSkipUpSolve_StateMachine
      end
    end

    # nope, regular clique init-solve, go to 1
    return isCliqUpSolved_StateMachine
end

"""
    $SIGNATURES

Notify possible parent if clique is upsolved and exit the state machine.

Notes
- State machine function nr.0
- can recycle if two checks:
  - previous clique was identically downsolved
  - all children are also :uprecycled
"""
function canCliqMargRecycle_StateMachine(csmc::CliqStateMachineContainer)
  # @show getCliqFrontalVarIds(csmc.oldcliqdata), getCliqueStatus(csmc.oldcliqdata)
  infocsm(csmc, "0., $(csmc.incremental) ? :uprecycled => getCliqueStatus(csmc.oldcliqdata)=$(getCliqueStatus(csmc.oldcliqdata))")

  if areCliqVariablesAllMarginalized(csmc.dfg, csmc.cliq)

    # no work required other than assembling upward message
    if getSolverParams(csmc.cliqSubFg).dbg
      tmnow = now()
      tmpst = getCliqueStatus(csmc.cliq)
      @async begin
        mkpath(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)"))
        open(joinLogPath(csmc.cliqSubFg,"logs","cliq$(csmc.cliq.index)","marginalization.log"), "w") do f
          println(f, tmnow, ", marginalized from previous status ", tmpst)
        end
      end
    end
    prepPutCliqueStatusMsgUp!(csmc, :marginalized, dfg=csmc.dfg)

    # set marginalized color
    setCliqDrawColor(csmc.cliq, "blue")

    # set flag, looks to be previously unused???
    getCliqueData(csmc.cliq).allmarginalized = true

    # FIXME divert to rapid CSM exit
	# GUESSING THIS THE RIGHT WAY go to 4
	# return canCliqMargSkipUpSolve_StateMachine
  end

  # go to 0c.
  return canCliqIncrRecycle_StateMachine
end







## =================================================================================================================


"""
    $SIGNATURES

Special internal function to try return the clique data if succesfully identified in `othertree::AbstractBayesTree`,
based on contents of `seeksSimilar::BayesTreeNodeData`.

Notes
- Used to identify and skip similar cliques (i.e. recycle computations)
"""
function attemptTreeSimilarClique(othertree::AbstractBayesTree, 
                                  seeksSimilar::BayesTreeNodeData  )
  #
  # inner convenience function for returning empty clique
  function EMPTYCLIQ()
    clq = TreeClique(-1,"null")
    setLabel!(clq, "")
    setCliqueData!(clq, BayesTreeNodeData())
    return clq
  end

  # does the other clique even exist?
  seekFrontals = getCliqFrontalVarIds(seeksSimilar)
  if !hasClique(othertree, seekFrontals[1])
    return EMPTYCLIQ()
  end

  # do the cliques share the same frontals?
  otherCliq = getClique(othertree, seekFrontals[1])
  otherFrontals = getCliqFrontalVarIds(otherCliq)
  commonFrontals = intersect(seekFrontals, otherFrontals)
  if length(commonFrontals) != length(seekFrontals) || length(commonFrontals) != length(otherFrontals)
    return EMPTYCLIQ()
  end

  # do the cliques share the same separator variables?
  seekSeparator = getCliqSeparatorVarIds(seeksSimilar)
  otherSeparator = getCliqSeparatorVarIds(otherCliq)
  commonSep = intersect(seekSeparator, otherSeparator)
  if length(commonSep) != length(seekSeparator) || length(commonSep) != length(otherSeparator)
    return EMPTYCLIQ()
  end

  # do the cliques use the same factors (potentials)
  seekPotentials = getCliqFactorIds(seeksSimilar)
  otherFactors = getCliqFactorIds(otherCliq)
  commonFactors = intersect(seekPotentials, otherFactors)
  if length(commonFactors) != length(seekPotentials) || length(commonFactors) != length(otherFactors)
    return EMPTYCLIQ()
  end

  # lets assume they are the same
  return otherCliq::TreeClique
end
