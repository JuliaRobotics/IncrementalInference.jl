
"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 3
"""
function solveUp_ParametricStateMachine_Old(csmc::CliqStateMachineContainer)
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
  if result.f_converged || result.g_converged
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

# solve relatives ignoring any priors keeping `from` at ϵ
# if clique has priors : solve to get a prior on `from`
# send messages as factors or just the beliefs? for now factors
function solveUp_ParametricStateMachine(csmc::CliqStateMachineContainer)
  infocsm(csmc, "Par-3, Solving Up")

  setCliqueDrawColor!(csmc.cliq, "red")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  msgfcts = Symbol[]

  for (idx, upmsg) in getMessageBuffer(csmc.cliq).upRx #get cached messages taken from children saved in this clique
    child_factors = addMsgFactors_Parametric!(csmc.cliqSubFg, upmsg, UpwardPass)
    append!(msgfcts, getLabel.(child_factors)) # addMsgFactors_Parametric!
  end
  logCSM(csmc, "length mgsfcts=$(length(msgfcts))")
  infocsm(csmc, "length mgsfcts=$(length(msgfcts))")

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_beforeupsolve")

  subfg = csmc.cliqSubFg

  frontals = getCliqFrontalVarIds(csmc.cliq)
  separators = getCliqSeparatorVarIds(csmc.cliq)

  # if its a root do full solve
  if length(getParent(csmc.tree, csmc.cliq)) == 0
    # M, vartypeslist, lm_r, Σ = solve_RLM(subfg; is_sparse=false, finiteDiffCovariance=true)
    autoinitParametric!(subfg)
    M, vartypeslist, lm_r, Σ = solveGraphParametric!(subfg; is_sparse=false, finiteDiffCovariance=true, damping_term_min=1e-18)

  else

    # select first seperator as constant reference at the identity element
    isempty(separators) && @warn "empty separators solving cliq $(csmc.cliq.id.value)" ls(subfg) lsf(subfg)
    from = first(separators)
    from_v = getVariable(subfg, from)
    getSolverData(from_v, :parametric).val[1] = getPointIdentity(getVariableType(from_v))

    #TODO handle priors
    # Variables that are free to move
    free_vars = [frontals; separators[2:end]]
    # Solve for the free variables

    @assert !isempty(lsf(subfg)) "No factors in clique $(csmc.cliq.id.value) ls=$(ls(subfg)) lsf=$(lsf(subfg))"

    # M, vartypeslist, lm_r, Σ = solve_RLM_conditional(subfg, free_vars, [from];)
    M, vartypeslist, lm_r, Σ = solve_RLM_conditional(subfg, free_vars, [from]; finiteDiffCovariance=false, damping_term_min=1e-18)

  end
  
  # FIXME check solve convergence
  if !true
    @error "Par-3, clique $(csmc.cliq.id) failed to converge in upsolve" result
    # propagate error to cleanly exit all cliques
    putErrorUp(csmc)
    if length(getParent(csmc.tree, csmc.cliq)) == 0
      putErrorDown(csmc)
      return IncrementalInference.exitStateMachine
    end

    return waitForDown_StateMachine
  end

  logCSM(csmc, "$(csmc.cliq.id): subfg solve converged sending messages")
  
  # Pack results in massage factors
  
  sigmas = extractMarginalsAP(M, vartypeslist, Σ)
  
  # FIXME fix MsgRelativeType
  relative_message_factors = MsgRelativeType();
  for (i, to) in enumerate(vartypeslist)
    if to in separators
      #assume full dim factor
      factype = selectFactorType(subfg, from, to)
      # make S symetrical
      # S = sigmas[i] # FIXME for some reason SMatrix is not invertable even though it is!!!!!!!!
      S = Matrix(sigmas[i])# FIXME
      S = (S + S') / 2
      # @assert all(isapprox.(S, sigmas[i], rtol=1e-3)) "Bad covariance matrix - not symetrical"
      !all(isapprox.(S, sigmas[i], rtol=1e-3)) && @error("Bad covariance matrix - not symetrical")
      # @assert all(diag(S) .> 0) "Bad covariance matrix - not positive diag"
      !all(diag(S) .> 0) && @error("Bad covariance matrix - not positive diag")

      
      M_to = getManifold(getVariableType(subfg, to))
      ϵ = getPointIdentity(M_to)
      μ = vee(M_to, ϵ, log(M_to, ϵ, lm_r[i]))
        
      message_factor = AdFactor(factype(MvNormal(μ, S)))
      
     
      # logCSM(csmc, "$(csmc.cliq.id): Z=$(getMeasurementParametric(message_factor))"; loglevel = Logging.Warn)

      push!(relative_message_factors, (variables=[from, to], likelihood=message_factor))
    end
  end

  # Done with solve delete factors
  #TODO confirm, maybe don't delete mesage factors on subgraph, maybe delete if its priors, but not conditionals
  # deleteMsgFactors!(csmc.cliqSubFg)

  # store the cliqSubFg for later debugging
  _dbgCSMSaveSubFG(csmc, "fg_afterupsolve")

  # cliqueLikelihood = calculateMarginalCliqueLikelihood(vardict, Σ, varIds, cliqSeparatorVarIds)

  #Fill in CliqueLikelihood
  beliefMsg = LikelihoodMessage(;
    sender = (; id = csmc.cliq.id.value, step = csmc._csm_iter),
    status = UPSOLVED,
    variableOrder = separators,
    # cliqueLikelihood,
    jointmsg = _MsgJointLikelihood(;relatives=relative_message_factors),
    msgType = ParametricMessage(),
  )

  # @assert length(separators) <= 2 "TODO length(separators) = $(length(separators)) > 2 in clique $(csmc.cliq.id.value)"
  @assert isempty(lsfPriors(csmc.cliqSubFg)) || csmc.cliq.id.value == 1 "TODO priors in clique $(csmc.cliq.id.value)"
  # if length(lsfPriors(csmc.cliqSubFg)) > 0 || length(separators) > 2
    # for si in cliqSeparatorVarIds
    #   vnd = getSolverData(getVariable(csmc.cliqSubFg, si), :parametric)
    #   beliefMsg.belief[si] = TreeBelief(deepcopy(vnd))
    # end
  # end

  for e in getEdgesParent(csmc.tree, csmc.cliq)
    logCSM(csmc, "$(csmc.cliq.id): put! on edge $(e)")
    getMessageBuffer(csmc.cliq).upTx = deepcopy(beliefMsg)
    putBeliefMessageUp!(csmc.tree, e, beliefMsg)
  end

  return waitForDown_StateMachine
end

global g_n = nothing

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
        
        p = belief.val[1]

        S = belief.bw
        S = (S + S') / 2
        vnd.bw .= S
        
        nd = MvNormal(getCoordinates(Main.Pose2, p), S)
        addFactor!(csmc.cliqSubFg, [msym], Main.PriorPose2(nd))
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
    # vardict, result, flatvars, Σ = solveConditionalsParametric(csmc.cliqSubFg, frontals)
    #TEMP testing difference
    # vardict, result = solveGraphParametric(csmc.cliqSubFg)
    # Pack all results in variables
    @assert !isempty(lsf(csmc.cliqSubFg)) "No factors in clique $(csmc.cliq.id.value) ls=$(ls(csmc.cliqSubFg)) lsf=$(lsf(csmc.cliqSubFg))"

    # M, vartypeslist, lm_r, Σ = solve_RLM_conditional(csmc.cliqSubFg, frontals; finiteDiffCovariance=false, damping_term_min=1e-18)
    M, vartypeslist, lm_r, Σ = solve_RLM(csmc.cliqSubFg; finiteDiffCovariance=false, damping_term_min=1e-18)
    sigmas = extractMarginalsAP(M, vartypeslist, Σ)

    if true # TODO check for convergence result.g_converged || result.f_converged
      logCSM(
        csmc,
        "$(csmc.cliq.id): subfg optim converged updating variables";
        loglevel = Logging.Debug,
      )
      for (i, v) in enumerate(vartypeslist)
        if v in frontals
          # logCSM(csmc, "$(csmc.cliq.id) down: updating $v"; val, loglevel = Logging.Debug)
          vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
          
          S = Matrix(sigmas[i])# FIXME
          S = (S + S') / 2
          # @assert all(isapprox.(S, sigmas[i], rtol=1e-3)) "Bad covariance matrix - not symetrical"
          !all(isapprox.(S, sigmas[i], rtol=1e-3)) && @error("Bad covariance matrix - not symetrical")
          # @assert all(diag(S) .> 0) "Bad covariance matrix - not positive diag"
          !all(diag(S) .> 0) && @error("Bad covariance matrix - not positive diag")

          
          #Update subfg variables
          vnd.val[1] = lm_r[i]
          vnd.bw .= S
        end
      end
      # for (v, val) in vardict
      #   logCSM(csmc, "$(csmc.cliq.id) down: updating $v"; val, loglevel = Logging.Debug)
      #   vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
        
      #   #Update subfg variables
      #   vnd.val[1] = val.val
      #   vnd.bw .= val.cov
      # end
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
    logCSM(csmc, "$(csmc.cliq.id): down message $fi"; beliefMsg=beliefMsg.belief[fi], loglevel = Logging.Debug)
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
