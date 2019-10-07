# clique state machine for tree based initialization and inference

"""
    $SIGNATURES

Specialized info logger print function to show clique state machine information
in a standardized form.
"""
function infocsm(csmc::CliqStateMachineContainer, str::A) where {A <: AbstractString}

  tm = string(Dates.now())
  tmt = split(tm, 'T')[end]

  lbl = csmc.cliq.attributes["label"]
  lbl1 = split(lbl,',')[1]
  cliqst = getCliqStatus(csmc.cliq)

  with_logger(csmc.logger) do
    @info "$tmt | $(csmc.cliq.index)---$lbl1 @ $(cliqst) | "*str
  end
  flush(csmc.logger.stream)
  nothing
end

"""
    $SIGNATURES

Do cliq downward inference

Notes:
- State machine function nr. 11
"""
function doCliqDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine")
  setCliqDrawColor(csmc.cliq, "red")
  opts = getSolverParams(csmc.dfg)
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(opts.logpath,"bt.pdf")) : nothing

  # get down msg from parent (assing root clique CSM wont make it here)
  prnt = getParent(csmc.tree, csmc.cliq)
  dwnmsgs = getDwnMsgs(prnt[1])
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- dwnmsgs=$(collect(keys(dwnmsgs)))")

  # maybe cycle through separators (or better yet, just use values directly -- see next line)
  msgfcts = addMsgFactors!(csmc.cliqSubFg, dwnmsgs)
  # force separator variables to adopt down message values
  updateSubFgFromDownMsgs!(csmc.cliqSubFg, dwnmsgs, getCliqSeparatorVarIds(csmc.cliq))

  # add required all frontal connected factors
  newvars, newfcts = addDownVariableFactors!(csmc.dfg, csmc.cliqSubFg, csmc.cliq, csmc.logger)

  # store the cliqSubFg for later debugging
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedownsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedownsolve.pdf"))
  end

  ## new way
  # calculate belief on each of the frontal variables and iterate if required
  solveCliqDownFrontalProducts!(csmc.cliqSubFg, csmc.cliq, opts, csmc.logger)

  # compute new down messages
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- going to set new down msgs.")
  getSetDownMessagesComplete!(csmc.cliqSubFg, csmc.cliq, dwnmsgs, csmc.logger)
  # setDwnMsg!(cliq, drt.keepdwnmsgs)

      # frsyms = getCliqFrontalVarIds(csmc.cliq)
      # # use new localproduct approach
      # # TODO refactor with @spawn jl 1.3
      # if csmc.opts.multiproc
      #   downresult = Dict{Symbol, Tuple{BallTreeDensity, Float64, Vector{Symbol}}}()
      #   @sync for i in 1:length(frsyms)
      #     @async begin
      #       downresult[frsyms[i]] = remotecall_fetch(localProductAndUpdate!, upp2(), csmc.cliqSubFg, frsyms[i])
      #     end
      #   end
      #   infocsm(csmc, "doCliqDownSolve_StateMachine, multiproc keys=$(keys(downresult))")
      #   for fr in frsyms
      #     infocsm(csmc, "doCliqDownSolve_StateMachine, multiproc, setValKDE! $fr, inferdim=$(downresult[fr][2]), lbls=$(downresult[fr][3])")
      #     setValKDE!(csmc.cliqSubFg, fr, downresult[fr][1], false, downresult[fr][2])
      #   end
      # else
      #   for fr in frsyms
      #     localProductAndUpdate!(csmc.cliqSubFg, fr, csmc.logger)
      #   end
      # end

    ## old approach
      # # multiproc = false
      # # call down inference, TODO multiproc
      # if csmc.opts.multiproc
      #   cliqc = deepcopy(csmc.cliq)
      #   cliqcd = getData(cliqc)
      #   # redirect to new unused so that CAN be serialized
      #   cliqcd.initUpChannel = Channel{Symbol}(1)
      #   cliqcd.initDownChannel = Channel{Symbol}(1)
      #   cliqcd.solveCondition = Condition()
      #   cliqcd.statehistory = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}()
      #   infocsm(csmc, "11, doCliqDownSolve_StateMachine -- MULTIPROC downGibbsCliqueDensity")
      #   # NOTE Cannot send logger to separate process memory
      #   drt = remotecall_fetch(downGibbsCliqueDensity, upp2(), csmc.cliqSubFg, cliqc, dwnmsgs, 100, 3, false, true)
      # else
      #   infocsm(csmc, "11, doCliqDownSolve_StateMachine -- SINGLEPROC downGibbsCliqueDensity")
      #   # usemsgpriors = true so that old plumbing is used properly
      #   drt = downGibbsCliqueDensity(csmc.cliqSubFg, csmc.cliq, dwnmsgs, 100, 3, false, true, csmc.logger)
      # end

      # # update clique subgraph with new status
      # updateFGBT!(csmc.cliqSubFg, csmc.tree, csmc.cliq.index, drt, dbg=false, fillcolor="lightblue", logger=csmc.logger)
      setCliqDrawColor(csmc.cliq, "lightblue")

  csmc.dodownsolve = false
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- finished with downGibbsCliqueDensity, now update csmc")

  # set MAP est
  for sym in getCliqFrontalVarIds(csmc.cliq)
    setVariablePosteriorEstimates!(csmc.cliqSubFg, sym)
  end

  # store the cliqSubFg for later debugging
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterdownsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterdownsolve.pdf"))
  end

  # transfer results to main factor graph
  frsyms = getCliqFrontalVarIds(csmc.cliq)
  infocsm(csmc, "11, finishingCliq -- going for transferUpdateSubGraph! on $frsyms")
  transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger)

  # setCliqStatus!(csmc.cliq, :downsolved) # should be a notify
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- before notifyCliqDownInitStatus!")
  notifyCliqDownInitStatus!(csmc.cliq, :downsolved, logger=csmc.logger)
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- just notified notifyCliqDownInitStatus!")

  # remove msg factors that were added to the subfg
  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- removing up message factors, length=$(length(msgfcts))")
  deleteMsgFactors!(csmc.cliqSubFg, msgfcts)

  # show tree if required
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(opts.logpath,"bt.pdf")) : nothing

  infocsm(csmc, "11, doCliqDownSolve_StateMachine -- finished, exiting CSM on clique=$(csmc.cliq.index)")

  # and finished
  return IncrementalInference.exitStateMachine
end


"""
    $SIGNATURES

Direct state machine to continue with downward solve or exit.

Notes
- State machine function nr. 10
"""
function determineCliqIfDownSolve_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "10, determineCliqIfDownSolve_StateMachine, csmc.dodownsolve=$(csmc.dodownsolve).")
  # finished and exit downsolve
  if !csmc.dodownsolve
    infocsm(csmc, "10, determineCliqIfDownSolve_StateMachine -- shortcut exit since downsolve not required.")
    return IncrementalInference.exitStateMachine
  end

  # yes, continue with downsolve
  setCliqDrawColor(csmc.cliq, "turquoise")
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  # assume separate down solve via solveCliq! call, but need a csmc.cliqSubFg this late in CSM anyway -- so just go copy one
  if length(ls(csmc.cliqSubFg)) == 0
    # first need to fetch cliq sub graph
    # go to 2b
    return buildCliqSubgraphForDown_StateMachine
  end

  # block here until parent is downsolved
  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    infocsm(csmc, "10, determineCliqIfDownSolve_StateMachine, going to block on parent.")
    # TODO -- some cleanup
    blockCliqUntilParentDownSolved(prnt[1], logger=csmc.logger)
    prntst = getCliqStatus(prnt[1])
    infocsm(csmc, "10, determineCliqIfDownSolve_StateMachine, parent status=$prntst.")
    if prntst != :downsolved
      infocsm(csmc, "10, determineCliqIfDownSolve_StateMachine, going around again.")
      return determineCliqIfDownSolve_StateMachine
  	end
  else
    # special case for down solve on root clique.  When using solveCliq! following an up pass.

    # this is the root clique, so assume already downsolved -- only special case
    dwnmsgs = getCliqDownMsgsAfterDownSolve(csmc.cliqSubFg, csmc.cliq)
    setCliqDrawColor(csmc.cliq, "lightblue")
    setDwnMsg!(csmc.cliq, dwnmsgs)
    setCliqStatus!(csmc.cliq, :downsolved)
	csmc.dodownsolve = false

	# Update estimates and transfer back to the graph
	frsyms = getCliqFrontalVarIds(csmc.cliq)
	# Calculate estimates
	map(sym -> setVariablePosteriorEstimates!(csmc.cliqSubFg, sym), frsyms)
	# Transfer to parent graph
	transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms)

    notifyCliqDownInitStatus!(csmc.cliq, :downsolved, logger=csmc.logger)

    return IncrementalInference.exitStateMachine
  end

  infocsm(csmc, "10, going for down solve.")
  # go to 11
  return doCliqDownSolve_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr.9
"""
function finishCliqSolveCheck_StateMachine(csmc::CliqStateMachineContainer)
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing
  cliqst = getCliqStatus(csmc.cliq)
  infocsm(csmc, "9, finishingCliq")
  if cliqst == :upsolved
      frsyms = getCliqFrontalVarIds(csmc.cliq)
    infocsm(csmc, "9, finishingCliq -- going for transferUpdateSubGraph! on $frsyms")
    # TODO what about down solve??
    transferUpdateSubGraph!(csmc.dfg, csmc.cliqSubFg, frsyms, csmc.logger)

    # remove any solvable upward cached data -- TODO will have to be changed for long down partial chains
    # assuming maximally complte up solved cliq at this point
    lockUpStatus!(csmc.cliq)
    sdims = Dict{Symbol,Float64}()
    for varid in getCliqAllVarIds(csmc.cliq)
      sdims[varid] = 0.0
    end
    updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
    unlockUpStatus!(csmc.cliq)

    # go to 10
    return determineCliqIfDownSolve_StateMachine # IncrementalInference.exitStateMachine
  elseif cliqst == :initialized
    setCliqDrawColor(csmc.cliq, "sienna")
    csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  else
    infocsm(csmc, "9, finishingCliq -- init not complete and should wait on init down message.")
    setCliqDrawColor(csmc.cliq, "coral")
    # TODO, potential problem with trying to downsolve
    # return isCliqNull_StateMachine # doesCliqNeeddownmsg_StateMachine
  end

  # go to 4
  return isCliqNull_StateMachine # whileCliqNotSolved_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr. 8c
"""
function waitChangeOnParentCondition_StateMachine(csmc::CliqStateMachineContainer)
  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, wait on parent=$(prnt[1].index) for condition notify.")
    @sync begin
      @async begin
        sleep(1)
        notify(getSolveCondition(prnt[1]))
      end
      wait(getSolveCondition(prnt[1]))
    end
  else
    infocsm(csmc, "8c, waitChangeOnParentCondition_StateMachine, cannot wait on parent for condition notify.")
    @warn "no parent!"
  end

  # go to 4
  return isCliqNull_StateMachine
end


"""
    $SIGNATURES

Do up initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in upward direction.

Notes
- State machine function nr. 8b
- Includes initialization routines.
- TODO: Make multi-core
"""
function attemptCliqInitUp_StateMachine(csmc::CliqStateMachineContainer)
  # get a few pointers to relevant data
  cliqst = getCliqStatus(csmc.cliq)
  opts = getSolverParams(csmc.dfg)

  infocsm(csmc, "8b, attemptCliqInitUp, !areCliqChildrenNeedDownMsg()=$(!areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq))" )
  if cliqst in [:initialized; :null; :needdownmsg] && !areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
    setCliqDrawColor(csmc.cliq, "red")
    csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(opts.logpath,"bt.pdf")) : nothing

    # check if init is required and possible
    infocsm(csmc, "8b, attemptCliqInitUp, going for doCliqAutoInitUpPart1!.")
    # get incoming clique up messages
    upmsgs = getCliqInitUpMsgs(csmc.cliq)
    # add incoming up messages as priors to subfg
    infocsm(csmc, "8b, doCliqAutoInitUpPart1! -- adding up message factors")
    msgfcts = addMsgFactors!(csmc.cliqSubFg, upmsgs)

    # store the cliqSubFg for later debugging
    if opts.dbg
      DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve"))
      drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve.pdf"))
    end

    doCliqAutoInitUpPart1!(csmc.cliqSubFg, csmc.tree, csmc.cliq, logger=csmc.logger)
    infocsm(csmc, "8b, attemptCliqInitUp, areCliqVariablesAllInitialized(subfg, cliq)=$(areCliqVariablesAllInitialized(csmc.cliqSubFg, csmc.cliq))")

    # do actual up solve
    retstatus = doCliqAutoInitUpPart2!(csmc.cliqSubFg, csmc.tree, csmc.cliq, multiproc=csmc.opts.multiproc, logger=csmc.logger)

    # remove msg factors that were added to the subfg
    infocsm(csmc, "8b, doCliqAutoInitUpPart2! -- removing up message factors, length=$(length(msgfcts))")
    deleteMsgFactors!(csmc.cliqSubFg, msgfcts)

    # store the cliqSubFg for later debugging
    if opts.dbg
      DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve"))
      drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve.pdf"))
    end

    # notify of results
    if cliqst != retstatus
      infocsm(csmc, "8b, attemptCliqInitUp, post-doCliqAu. -- notification retstatus=$retstatus")
      notifyCliqUpInitStatus!(csmc.cliq, retstatus, logger=csmc.logger)
    else
      infocsm(csmc, "8b, attemptCliqInitUp, post-doCliqAu. -- no notification required $cliqst=$retstatus")
    end
  end

  # go to 9
  return finishCliqSolveCheck_StateMachine
end


"""
    $SIGNATURES

Do down initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in downward direction.

Notes
- State machine function nr. 8a
- Includes initialization routines.
- TODO: Make multi-core
"""
function attemptCliqInitDown_StateMachine(csmc::CliqStateMachineContainer)
  #
  # should never happen to
  setCliqDrawColor(csmc.cliq, "green")
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  # initialize clique in downward direction
  # not if parent also needs downward init message
  infocsm(csmc, "8a, needs down message -- attempt down init")
  prnt = getParent(csmc.tree, csmc.cliq)[1]

  # take atomic lock when waiting for down ward information
  lockUpStatus!(getData(prnt))

  dwinmsgs = prepCliqInitMsgsDown!(csmc.dfg, csmc.tree, prnt, csmc.cliq, logger=csmc.logger, dbgnew=!haskey(getSolverParams(csmc.dfg).devParams,:dontUseParentFactorsInitDown)) # csmc.cliqSubFg
  dwnkeys = collect(keys(dwinmsgs))


  infocsm(csmc, "8a, attemptCliqInitD., dwinmsgs=$(dwnkeys), adding msg factors")

  ## DEVIdea
  msgfcts = addMsgFactors!(csmc.cliqSubFg, dwinmsgs)

  # determine if more info is needed for partial

  # DEVidea
  sdims = getCliqVariableMoreInitDims(csmc.cliqSubFg, csmc.cliq)
  updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)


  infocsm(csmc, "8a, attemptCliqInitD., updated clique solvable dims")

  # priorize solve order for mustinitdown with lowest dependency first
  # follow example from issue #344
  mustwait = false
  if length(intersect(dwnkeys, getCliqSeparatorVarIds(csmc.cliq))) == 0 # length(dwinmsgs) == 0 ||
    infocsm(csmc, "8a, attemptCliqInitDown_StateMachine, no can do, must wait for siblings to update parent first.")
    mustwait = true
  elseif getSiblingsDelayOrder(csmc.tree, csmc.cliq, prnt, dwinmsgs, logger=csmc.logger)
    infocsm(csmc, "8a, attemptCliqInitD., prioritize")
    mustwait = true
  elseif getCliqSiblingsPartialNeeds(csmc.tree, csmc.cliq, prnt, dwinmsgs, logger=csmc.logger)
    infocsm(csmc, "8a, attemptCliqInitD., partialneedsmore")
    mustwait = true
  end

  # remove the downward messages too
  deleteMsgFactors!(csmc.cliqSubFg, msgfcts)

  infocsm(csmc, "8a, attemptCliqInitD., deleted msg factors and unlockUpStatus!")
  # unlock
  unlockUpStatus!(getData(prnt))
  infocsm(csmc, "8a, attemptCliqInitD., unlocked")

  solord = getCliqSiblingsPriorityInitOrder( csmc.tree, prnt, csmc.logger )
  infocsm(csmc, "8a, attemptCliqInitD., before are Siblings")
  noOneElse = areSiblingsRemaingNeedDownOnly(csmc.tree, csmc.cliq)
  infocsm(csmc, "8a, attemptCliqInitD., $mustwait, $noOneElse, solord = $solord")

  if mustwait && csmc.cliq.index!=solord[1] # && !noOneElse
    # go to 8c
    infocsm(csmc, "8a, attemptCliqInitD., mustwait, so wait on change.")
    return waitChangeOnParentCondition_StateMachine
  end

  ## TODO deal with partial inits only, either delay or continue at end...
  # find intersect between downinitmsgs and local clique variables
  # if only partials available, then

  infocsm(csmc, "8a, attemptCliqInitD.,do cliq init down dwinmsgs=$(keys(dwinmsgs))")
  cliqst = doCliqInitDown!(csmc.cliqSubFg, csmc.cliq, dwinmsgs, dbg=getSolverParams(csmc.dfg).dbg, logger=csmc.logger, logpath=getSolverParams(csmc.dfg).logpath )
  # TODO: transfer values changed in the cliques should be transfered to the tree in proc 1 here.

  # # TODO: is status of notify required here?
  setCliqStatus!(csmc.cliq, cliqst)
  # notifyCliqUpInitStatus!(csmc.cliq, cliqst)

  # TODO move out
  children = getChildren(csmc.tree, csmc.cliq)
  if areCliqChildrenNeedDownMsg(children)
    # set messages if children :needdownmsg
    infocsm(csmc, "8a, doCliqInitDown! -- must set messages for future down init")
    # construct init's up msg to place in parent from initialized separator variables
    msg = prepCliqInitMsgsUp(csmc.cliqSubFg, csmc.cliq, csmc.logger) # , tree,

    infocsm(csmc, "8a, putting fake upinitmsg in this cliq, msgs labels $(collect(keys(msg)))")
    # set fake up and notify down status -- repeat change status to same as notifyUp above
    setCliqUpInitMsgs!(csmc.cliq, csmc.cliq.index, msg)
    # setCliqStatus!(csmc.cliq, cliqst)
    setCliqDrawColor(csmc.cliq, "sienna")
    csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

    notifyCliqDownInitStatus!(csmc.cliq, cliqst, logger=csmc.logger)

    infocsm(csmc, "8a, near-end down init attempt, $cliqst.")
  end

  return attemptCliqInitUp_StateMachine
end


"""
    $SIGNATURES

Delay loop if waiting on upsolves to complete.

Notes
- State machine 7b
"""
function slowCliqIfChildrenNotUpsolved_StateMachine(csmc::CliqStateMachineContainer)

  # special case short cut
  cliqst = getCliqStatus(csmc.cliq)
  if cliqst == :needdownmsg
    infocsm(csmc, "7b, shortcut on cliq is a needdownmsg status.")
    return isCliqNull_StateMachine
  end
  childs = getChildren(csmc.tree, csmc.cliq)
  len = length(childs)
  @inbounds for i in 1:len
    if !(getCliqStatus(childs[i]) in [:upsolved;:uprecycled;:marginalized])
      infocsm(csmc, "7b, wait condition on upsolve, i=$i, ch_lbl=$(getCliqFrontalVarIds(childs[i])[1]).")
      wait(getSolveCondition(childs[i]))
      break
    end
  end

  # go to 4
  return isCliqNull_StateMachine
end

"""
    $SIGNATURES

Notes
- State machine function nr. 7
"""
function determineCliqNeedDownMsg_StateMachine(csmc::CliqStateMachineContainer)

  # fetch children status
  stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq, csmc.logger)

  # hard assumption here on upsolve from leaves to root
  proceed = true
  # fetch status from children (should already be available -- i.e. should not block)
  for (clid, clst) in stdict
    infocsm(csmc, "7, check stdict children: clid=$(clid), clst=$(clst)")
    # :needdownmsg # 'send' downward init msg direction
    !(clst in [:initialized;:upsolved;:marginalized;:downsolved;:uprecycled]) ? (proceed = false) : nothing
  end
  infocsm(csmc, "7, proceed=$(proceed), forceproceed=$(csmc.forceproceed)")


  if proceed || csmc.forceproceed
    csmc.forceproceed = false
    # return doCliqInferAttempt_StateMachine
    cliqst = getCliqStatus(csmc.cliq)
    infocsm(csmc, "7, status=$(cliqst), before attemptCliqInitDown_StateMachine")
    # d1,d2,cliqst = doCliqInitUpOrDown!(csmc.cliqSubFg, csmc.tree, csmc.cliq, isprntnddw)
    if cliqst == :needdownmsg && !isCliqParentNeedDownMsg(csmc.tree, csmc.cliq, csmc.logger)
      return attemptCliqInitDown_StateMachine
    end
    return attemptCliqInitUp_StateMachine
  else
    return slowCliqIfChildrenNotUpsolved_StateMachine
  end
end

"""
    $SIGNATURES

Notes
- State machine function nr. 6c
"""
function blockCliqSiblingsParentChildrenNeedDown_StateMachine(csmc::CliqStateMachineContainer)
  # add blocking case when all siblings and parent :needdownmsg -- until parent :initialized
  infocsm(csmc, "6c, check/block sibl&prnt :needdownmsg")
  blockCliqSiblingsParentNeedDown(csmc.tree, csmc.cliq, logger=csmc.logger)

  return determineCliqNeedDownMsg_StateMachine
end


"""
    $SIGNATURES

Notes
- State machine function nr. 5
"""
function blockUntilSiblingsStatus_StateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "5, blocking on parent until all sibling cliques have valid status")
  setCliqDrawColor(csmc.cliq, "blueviolet")
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  cliqst = getCliqStatus(csmc.cliq)
  infocsm(csmc, "5, block on siblings")
  prnt = getParent(csmc.tree, csmc.cliq)
  if length(prnt) > 0
    infocsm(csmc, "5, has parent clique=$(prnt[1].index)")
    blockCliqUntilChildrenHaveUpStatus(csmc.tree, prnt[1], csmc.logger)
  end

  # go to 6c
  return blockCliqSiblingsParentChildrenNeedDown_StateMachine
end



"""
    $SIGNATURES

Notes
- State machine function nr. 4
"""
function isCliqNull_StateMachine(csmc::CliqStateMachineContainer)

  prnt = getParent(csmc.tree, csmc.cliq)
  infocsm(csmc, "4, isCliqNull_StateMachine, csmc.incremental=$(csmc.incremental), len(prnt)=$(length(prnt))")
  #must happen before if :null
  stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq, csmc.logger)
  csmc.forceproceed = false

  # for recycle computed clique values case
  if csmc.incremental && getCliqStatus(csmc.oldcliqdata) == :downsolved
    csmc.incremental = false
    # might be able to recycle the previous clique solve, go to 0b
    return checkChildrenAllUpRecycled_StateMachine
  end

  if 0 == length(prnt)
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  end

  # go to 4b
  return doesCliqNeeddownmsg_StateMachine
end


"""
    $SIGNATURES

Determine if any down messages are required.

Notes
- State machine function nr.4b
"""
function doesCliqNeeddownmsg_StateMachine(csmc::CliqStateMachineContainer)

  cliqst = getCliqStatus(csmc.cliq)
  infocsm(csmc, "4b, doesCliqNeeddownmsg_StateMachine, cliqst=$cliqst")

  if cliqst != :null
    if cliqst != :needdownmsg
      return blockCliqSiblingsParentChildrenNeedDown_StateMachine
    end
  else
    # fetch (should not block)
    stdict = blockCliqUntilChildrenHaveUpStatus(csmc.tree, csmc.cliq, csmc.logger)
    chstatus = collect(values(stdict))
    len = length(chstatus)

    # if all children needdownmsg
    if len > 0 && sum(chstatus .== :needdownmsg) == len
      # TODO maybe can happen where some children need more information?
      infocsm(csmc, "4b, escalating to :needdownmsg since all children :needdownmsg")
      notifyCliqUpInitStatus!(csmc.cliq, :needdownmsg, logger=csmc.logger)
      setCliqDrawColor(csmc.cliq, "yellowgreen")
      csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

      # go to 5
      return blockUntilSiblingsStatus_StateMachine
    end

    # go to 6c
    return blockCliqSiblingsParentChildrenNeedDown_StateMachine
  end # != :null

  # if cliqst == :needdownmsg
    if areCliqChildrenNeedDownMsg(csmc.tree, csmc.cliq)
      infocsm(csmc, "4, must deal with child :needdownmsg")
      csmc.forceproceed = true
    else
      # go to 5
      return blockUntilSiblingsStatus_StateMachine
    end
  # end

  # go to 6c
  return blockCliqSiblingsParentChildrenNeedDown_StateMachine
end


"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- State machine function nr.2
"""
function buildCliqSubgraph_StateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarIds(csmc.cliq)
  # NOTE add all frontal factor neighbors DEV CASE -- use getData(cliq).dwnPotentials instead
  # fnsyms = getCliqVarsWithFrontalNeighbors(csmc.dfg, csmc.cliq)

  infocsm(csmc, "2, build subgraph syms=$(syms)")
  csmc.cliqSubFg = buildSubgraphFromLabels(csmc.dfg, syms)

  # store the cliqSubFg for later debugging
  opts = getSolverParams(csmc.dfg)
  if opts.dbg
    mkpath(joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)"))
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build.pdf"))
  end

  # go to 4
  return isCliqNull_StateMachine
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
  csmc.cliqSubFg = buildSubgraphFromLabels(csmc.dfg, syms)

  opts = getSolverParams(csmc.dfg)
  # store the cliqSubFg for later debugging
  if opts.dbg
    mkpath(joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)"))
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build_down"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build_down.pdf"))
  end

  # go to 10
  return determineCliqIfDownSolve_StateMachine
end

"""
    $SIGNATURES

Notify possible parent if clique is upsolved and exit the state machine.

Notes
- State machine function nr.1
"""
function isCliqUpSolved_StateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "1, isCliqUpSolved_StateMachine")
  cliqst = getCliqStatus(csmc.cliq)

  if cliqst in [:upsolved; :downsolved; :marginalized; :uprecycled]  #moved to 4 --- csmc.incremental &&
    # prep and send upward message
    prnt = getParent(csmc.tree, csmc.cliq)
    if length(prnt) > 0
      # not a root clique
      # construct init's up msg to place in parent from initialized separator variables
      msg = prepCliqInitMsgsUp(csmc.dfg, csmc.tree, csmc.cliq)
      setCliqUpInitMsgs!(prnt[1], csmc.cliq.index, msg)
      notifyCliqUpInitStatus!(csmc.cliq, cliqst, logger=csmc.logger)
    end
    #go to 10
    return determineCliqIfDownSolve_StateMachine
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
    chst = getCliqStatus(ch)
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
    updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)
    setCliqStatus!(csmc.cliq, :uprecycled)
    setCliqDrawColor(csmc.cliq, "orange")
    csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing
    # go to 1
    return isCliqUpSolved_StateMachine
  end

  # return to regular solve, go to 2
  return buildCliqSubgraph_StateMachine
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
function testCliqCanRecycled_StateMachine(csmc::CliqStateMachineContainer)
  # @show getCliqFrontalVarIds(csmc.oldcliqdata), getCliqStatus(csmc.oldcliqdata)
  infocsm(csmc, "0., $(csmc.incremental) ? :uprecycled => getCliqStatus(csmc.oldcliqdata)=$(getCliqStatus(csmc.oldcliqdata))")

  # check if should be trying and can recycle clique computations
  if csmc.incremental && getCliqStatus(csmc.oldcliqdata) == :downsolved
    # check if a subgraph will be needed later
    if csmc.dodownsolve
      # yes need subgraph and need more checks, so go to 2
      return buildCliqSubgraph_StateMachine
    else
       # one or two checks say yes, so go to 4
      return isCliqNull_StateMachine
    end
  end

  # nope, regular clique init-solve, go to 1
  return isCliqUpSolved_StateMachine
end


"""
    $SIGNATURES

EXPERIMENTAL: perform upward inference using a state machine solution approach.

Notes:
- will call on values from children or parent cliques
- can be called multiple times
- Assumes all cliques in tree are being solved simultaneously and in similar manner.
- State machine rev.1 -- copied from first TreeBasedInitialization.jl.
- Doesn't do partial initialized state properly yet.
"""
function cliqInitSolveUpByStateMachine!(dfg::G,
                                        tree::BayesTree,
                                        cliq::Graphs.ExVertex;
                                        N::Int=100,
										oldcliqdata::BayesTreeNodeData=emptyBTNodeData(),
                                        drawtree::Bool=false,
                                        show::Bool=false,
                                        incremental::Bool=true,
                                        limititers::Int=-1,
                                        upsolve::Bool=true,
                                        downsolve::Bool=true,
                                        recordhistory::Bool=false,
                                        delay::Bool=false,
                                        logger::SimpleLogger=SimpleLogger(Base.stdout)) where {G <: AbstractDFG, AL <: AbstractLogger}
  #
  children = Graphs.ExVertex[]
  for ch in Graphs.out_neighbors(cliq, tree.bt)
    push!(children, ch)
  end
  prnt = getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : InMemDFGType#GraphsDFG{SolverParams}

  csmc = CliqStateMachineContainer(dfg, initfg(destType), tree, cliq, prnt, children, false, incremental, drawtree, downsolve, delay, getSolverParams(dfg), oldcliqdata, logger)

  nxt = upsolve ? testCliqCanRecycled_StateMachine : (downsolve ? testCliqCanRecycled_StateMachine : error("must attempt either up or down solve"))

  statemachine = StateMachine{CliqStateMachineContainer}(next=nxt)
  while statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory); end
  statemachine.history
end



#
