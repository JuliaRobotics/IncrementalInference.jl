
"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 3
"""
function solveUp_ParametricStateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "Par-3, Solving Up")

  setCliqueDrawColor!(csmc.cliq, "red")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  #TODO maybe change to symbols
  msgfcts = DFGFactor[]
  # LITTLE WEIRD get previously set up msgs (stored in this clique)
  # FIXME, fetch message buffered in channels
  # see #855
  for (idx, upmsg) in getMessageBuffer(csmc.cliq).upRx #get cached messages taken from children saved in this clique
    #TODO remove temp msgfcts container
    append!(msgfcts, addMsgFactors!(csmc.cliqSubFg, upmsg, UpwardPass)) # addMsgFactors_Parametric!
  end
  logCSM(csmc, "length mgsfcts=$(length(msgfcts))")
  infocsm(csmc, "length mgsfcts=$(length(msgfcts))")

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")

  vardict, result, varIds, Σ = solveGraphParametricOptim(csmc.cliqSubFg)

  logCSM(csmc, "$(csmc.cliq.id) vars $(keys(varIds))")
  # @info "$(csmc.cliq.id) Σ $(Σ)"
  # Pack all results in variables
  # FIXME test f_converged, ls_success, confirm convergence check
  if Optim.g_converged(result) || Optim.f_converged(result)
    logCSM(csmc, "$(csmc.cliq.id): subfg optim converged updating variables")
    for (v, val) in vardict
      vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
      # fill in the variable node data value
      logCSM(csmc, "$(csmc.cliq.id) up: updating $v : $val")
      vnd.val[1] = val.val
      #calculate and fill in covariance
      #TODO rather broadcast than make new memory
      vnd.bw = val.cov
    end
    # elseif length(lsfPriors(csmc.cliqSubFg)) == 0 #FIXME
    #   @error "Par-3, clique $(csmc.cliq.id) failed to converge in upsolve, but ignoring since no priors" result
  else
    @error "Par-3, clique $(csmc.cliq.id) failed to converge in upsolve" result
    # propagate error to cleanly exit all cliques
    putErrorUp(csmc)
    if length(getParent(csmc.tree, csmc.cliq)) == 0
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end

    return waitForDown_StateMachine
  end

  # Done with solve delete factors
  #TODO confirm, maybe don't delete mesage factors on subgraph, maybe delete if its priors, but not conditionals
  deleteMsgFactors!(csmc.cliqSubFg)

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterupsolve")

  #fill in belief
  #TODO createBeliefMessageParametric(csmc.cliqSubFg, csmc.cliq, solvekey=opts.solvekey)
  cliqSeparatorVarIds = getCliqSeparatorVarIds(csmc.cliq)
  #Fill in CliqueLikelihood
  cliqlikelihood =
    calculateMarginalCliqueLikelihood(vardict, Σ, varIds, cliqSeparatorVarIds)
  # @info "$(csmc.cliq.id) clique likelihood message $(cliqlikelihood)"
  beliefMsg = LikelihoodMessage(;
    sender = (; id = csmc.cliq.id.value, step = csmc._csm_iter),
    status = UPSOLVED,
    variableOrder = cliqSeparatorVarIds,
    cliqueLikelihood = cliqlikelihood,
    msgType = ParametricMessage(),
  )

  #FIXME bit of a hack, only fill in variable beliefs if there are priors or for now more than one seperator
  if length(lsfPriors(csmc.cliqSubFg)) > 0 || length(cliqSeparatorVarIds) > 1
    for si in cliqSeparatorVarIds
      vnd = getSolverData(getVariable(csmc.cliqSubFg, si), :parametric)
      beliefMsg.belief[si] = TreeBelief(deepcopy(vnd))
    end
  end

  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.id): put! on edge $(e)")
    getMessageBuffer(csmc.cliq).upTx = deepcopy(beliefMsg)
    putBeliefMessageUp!(csmc.tree, e, beliefMsg)
  end

  return waitForDown_StateMachine
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 5
"""
function solveDown_ParametricStateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "Par-5, Solving down")

  setCliqueDrawColor!(csmc.cliq, "red")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  # TODO create function: 
  # updateMsgSeparators!(csmc.cliqSubFg, downmsg)
  downmsg = getMessageBuffer(csmc.cliq).downRx  #see #855
  svars = getCliqSeparatorVarIds(csmc.cliq)
  if !isnothing(downmsg)
    for (msym, belief) in downmsg.belief
      if msym in svars
        #TODO maybe combine variable and factor in new prior?
        vnd = getSolverData(getVariable(csmc.cliqSubFg, msym), :parametric)
        logCSM(csmc, "$(csmc.cliq.id): Updating separator $msym from message $(belief.val)")
        vnd.val .= belief.val
        vnd.bw .= belief.bw
      end
    end
  end

  # store the cliqSubFg for later debugging
  # NOTE ITS not changed for now but keep here for possible future use
  # _dbgCSMSaveSubFG(csmc, "fg_beforedownsolve")

  # DownSolve cliqSubFg
  #only down solve if its not a root
  if length(getParent(csmc.tree, csmc.cliq)) != 0
    frontals = getCliqFrontalVarIds(csmc.cliq)
    vardict, result, flatvars, Σ = solveConditionalsParametric(csmc.cliqSubFg, frontals)
    #TEMP testing difference
    # vardict, result = solveGraphParametric(csmc.cliqSubFg)
    # Pack all results in variables
    if Optim.g_converged(result) || Optim.f_converged(result)
      logCSM(
        csmc,
        "$(csmc.cliq.id): subfg optim converged updating variables";
        loglevel = Logging.Info,
      )
      for (v, val) in vardict
        logCSM(csmc, "$(csmc.cliq.id) down: updating $v : $val"; loglevel = Logging.Info)
        vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
        #Update subfg variables
        vnd.val[1] = val.val
        vnd.bw .= val.cov
      end
    else
      @error "Par-5, clique $(csmc.cliq.id) failed to converge in down solve" result
      #propagate error to cleanly exit all cliques
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end
  end

  #TODO fill in belief
  cliqFrontalVarIds = getCliqFrontalVarIds(csmc.cliq)
  #TODO createBeliefMessageParametric
  # beliefMsg = createBeliefMessageParametric(csmc.cliqSubFg, cliqFrontalVarIds, solvekey=opts.solvekey)
  beliefMsg = LikelihoodMessage(;
    sender = (; id = csmc.cliq.id.value, step = csmc._csm_iter),
    status = DOWNSOLVED,
    msgType = ParametricMessage(),
  )
  for fi in cliqFrontalVarIds
    vnd = getSolverData(getVariable(csmc.cliqSubFg, fi), :parametric)
    beliefMsg.belief[fi] = TreeBelief(vnd)
    logCSM(csmc, "$(csmc.cliq.id): down message $fi : $beliefMsg"; loglevel = Logging.Info)
  end

  # pass through the frontal variables that were sent from above
  if !isnothing(downmsg)
    pass_through_separators = intersect(svars, keys(downmsg.belief))
    for si in pass_through_separators
      beliefMsg.belief[si] = downmsg.belief[si]
      logCSM(csmc, "adding parent message"; sym = si, msg = downmsg.belief[si])
    end
  end

  #TODO sendBeliefMessageParametric(csmc, beliefMsg)
  #TODO maybe send a specific message to only the child that needs it
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.id): put! on edge $(e)")
    @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messageChannels[e.index].downMsg, beliefMsg)
  end

  logCSM(csmc, "$(csmc.cliq.id): Solve completed")

  return updateFromSubgraph_StateMachine
end

#
