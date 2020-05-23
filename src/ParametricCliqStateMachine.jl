

"""
    $SIGNATURES

EXPERIMENTAL: Init and start state machine for parametric solve.
"""
function initStartCliqStateMachineParametric!(dfg::G,
                                              tree::AbstractBayesTree,
                                              cliq::TreeClique,
                                              cliqKey::Int;
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

  # TODO Taking out, use tree and messages for operations involving children and parents
  children = TreeClique[]#getChildren(tree, cliq)
  prnt = TreeClique[]#getParent(tree, cliq)

  destType = (G <: InMemoryDFGTypes) ? G : InMemDFGType#GraphsDFG{SolverParams}

  #csmc = CliqStateMachineContainer(dfg, initfg(destType), tree, cliq, prnt, children, false, incremental, drawtree, downsolve, delay, getSolverParams(dfg), oldcliqdata, logger)
  csmc = CliqStateMachineContainer(dfg, initfg(destType, solverParams=getSolverParams(dfg)), tree, cliq,
                                    cliqKey, prnt, children,
                                    false, incremental, drawtree, downsolve, delay,
                                    getSolverParams(dfg), Dict{Symbol,String}(), oldcliqdata, logger,
                                    Dict{Int,LikelihoodMessage}(), LikelihoodMessage())

  # nxt = upsolve ? testCliqCanRecycled_ParametricStateMachine : (downsolve ? testCliqCanRecycled_ParametricStateMachine : error("must attempt either up or down solve"))
  nxt = buildCliqSubgraph_ParametricStateMachine

  statemachine = StateMachine{CliqStateMachineContainer}(next=nxt)
  while statemachine(csmc, verbose=true, iterlimit=limititers, recordhistory=recordhistory); end
  statemachine.history
end

"""
    $SIGNATURES

Build a sub factor graph for clique variables from the larger factor graph.

Notes
- Parametric State machine function nr.1
"""
function buildCliqSubgraph_ParametricStateMachine(csmc::CliqStateMachineContainer)
  # build a local subgraph for inference operations
  syms = getCliqAllVarIds(csmc.cliq)

  infocsm(csmc, "Par-1, build subgraph syms=$(syms)")

  frontsyms = getCliqFrontalVarIds(csmc.cliq)
  sepsyms = getCliqSeparatorVarIds(csmc.cliq)
  buildCliqSubgraph!(csmc.cliqSubFg, csmc.dfg, frontsyms, sepsyms)

  #TODO remove, just as a sanity check if any priors remains on seperators
  removedIds = removeSeparatorPriorsFromSubgraph!(csmc.cliqSubFg, csmc.cliq)
  length(removedIds) > 0 && @error "removeSeparatorPriorsFromSubgraph, removed priors should not happen"


  # store the cliqSubFg for later debugging
  opts = getSolverParams(csmc.dfg)
  if opts.dbg
    mkpath(joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)"))
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_build.pdf"))
  end

  # go to 2 wait for up
  return waitForUp_ParametricStateMachine
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 2
"""
function waitForUp_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-2, wait for up messages of needed")

  setCliqDrawColor(csmc.cliq, "purple") #TODO don't know if this is correct color
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  childrenOk = true
  beliefMessages = Dict{Int,LikelihoodMessage}()
  # beliefMessages = LikelihoodMessage[]

  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    @async begin
      @info "$(csmc.cliq.index): take! on edge $(isa(e,Graphs.Edge) ? e.index : e)"
      # Blocks until data is available. -- pull model #674
      beliefMsg = takeBeliefMessageUp!(csmc.tree, e)#take!(csmc.tree.messages[e.index].upMsg)
      beliefMessages[e.target.index] = beliefMsg
      # push!(beliefMessages, beliefMsg)
      @info "$(csmc.cliq.index): Belief message received with status $(beliefMsg.status)"
    end
  end

  for (idx,beliefMsg) in beliefMessages
    #save up message (and add priors to cliqSubFg)
    #kies csmc vir boodskappe vir debugging, dis 'n vector een per kind knoop
    if beliefMsg.status == UPSOLVED
      setUpMsg!(csmc, idx, beliefMsg)
      # csmc.msgsUp[idx] = beliefMsg
      # push!(csmc.msgsUp, beliefMsg)

    else
      setCliqDrawColor(csmc.cliq, "red")
      # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

      for e in getEdgesParent(csmc.tree, csmc.cliq)
        @info "Par-2, $(csmc.cliq.index): propagate up error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        putBeliefMessageUp!(csmc.tree, e, LikelihoodMessage(ERROR_STATUS))#put!(csmc.tree.messages[e.index].upMsg,  LikelihoodMessage(ERROR_STATUS))
      end
      #if its the root, propagate error down
      #FIXME rather check if no parents with function (hasParents or isRoot)
      if length(getParent(csmc.tree, csmc.cliq)) == 0
        @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
          @info "Par-2 Root $(csmc.cliq.index): propagate down error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
          @async putBeliefMessageDown!(csmc.tree, e, LikelihoodMessage(ERROR_STATUS))#put!(csmc.tree.messages[e.index].downMsg,  LikelihoodMessage(ERROR_STATUS))
        end
        @error "Par-2 $(csmc.cliq.index): Exit with error state"
        return IncrementalInference.exitStateMachine
      end
      childrenOk = false
    end

  end

  !childrenOk && (return waitForDown_ParametricStateMachine)

  return solveUp_ParametricStateMachine
end

# Graph.jl does not have an in_edges function for a GenericIncidenceList, so extending here.
function Graphs.in_edges(vert::V, gr::GenericIncidenceList{V, Edge{V}, Vector{V}}) where {V}
  inclist = gr.inclist
  targid = vert.index
  inlist = Edge{V}[]
  for edgelist in inclist
    for ed in edgelist
      if ed.target.index == targid
        push!(inlist, ed)
      end
    end
  end
  return inlist
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 3
"""
function solveUp_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-3, Solving Up")

  setCliqDrawColor(csmc.cliq, "red")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  #TODO maybe change to symbols
  msgfcts = DFGFactor[]
  for (idx,upmsgs) in getUpMsgs(csmc) # csmc.msgsUp
    append!(msgfcts, addMsgFactors_Parametric!(csmc.cliqSubFg, upmsgs))
  end

  # store the cliqSubFg for later debugging
  opts = getSolverParams(csmc.dfg)
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve.pdf"))
  end

  vardict, result, varIds, Σ = solveFactorGraphParametric(csmc.cliqSubFg)

  @info "$(csmc.cliq.index) vars $(keys(varIds))"
  # @info "$(csmc.cliq.index) Σ $(Σ)"
  # Pack all results in variables
  # FIXME test f_converged, ls_success, confirm convergence check
  if result.f_converged || result.g_converged
    @info "$(csmc.cliq.index): subfg optim converged updating variables"
    for (v,val) in vardict
      vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
      # fill in the variable node data value
      @info "$(csmc.cliq.index) up: updating $v : $val"
      vnd.val .= val.val
      #calculate and fill in covariance
      #TODO rather broadcast than make new memory
      vnd.bw = val.cov
    end
  # elseif length(lsfPriors(csmc.cliqSubFg)) == 0 #FIXME
  #   @error "Par-3, clique $(csmc.cliq.index) failed to converge in upsolve, but ignoring since no priors" result
  else
    @error "Par-3, clique $(csmc.cliq.index) failed to converge in upsolve" result

    # propagate error to cleanly exit all cliques?
    beliefMsg = LikelihoodMessage(ERROR_STATUS)
    for e in getEdgesParent(csmc.tree, csmc.cliq)
      @info "Par-3, $(csmc.cliq.index): generate up error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
      putBeliefMessageUp!(csmc.tree, e, beliefMsg)# put!(csmc.tree.messages[e.index].upMsg, beliefMsg)
    end

    if length(getParent(csmc.tree, csmc.cliq)) == 0
      @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
        @info "Par-3 Root $(csmc.cliq.index): generate down error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        @async putBeliefMessageDown!(csmc.tree, e, LikelihoodMessage(ERROR_STATUS))#put!(csmc.tree.messages[e.index].downMsg,  LikelihoodMessage(ERROR_STATUS))
      end
      @error "Par-3 $(csmc.cliq.index): Exit with error state"
      return IncrementalInference.exitStateMachine
    end

    return waitForDown_ParametricStateMachine
  end

  # Done with solve delete factors
  #TODO confirm, maybe don't delete mesage factors on subgraph, maybe delete if its priors, but not conditionals
  deleteMsgFactors!(csmc.cliqSubFg, msgfcts)

  # store the cliqSubFg for later debugging
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve.pdf"))
  end

  #fill in belief
  #TODO createBeliefMessageParametric(csmc.cliqSubFg, csmc.cliq, solvekey=opts.solvekey)
  cliqSeparatorVarIds = getCliqSeparatorVarIds(csmc.cliq)
  #Fil in CliqueLikelihood
  cliqlikelihood = calculateMarginalCliqueLikelihood(vardict, Σ, varIds, cliqSeparatorVarIds)
  # @info "$(csmc.cliq.index) clique likelihood message $(cliqlikelihood)"
  beliefMsg = LikelihoodMessage(UPSOLVED, cliqSeparatorVarIds, cliqlikelihood)

  #FIXME bit of a hack, only fill in variable beliefs if there are priors or for now more than one seperator
  if length(lsfPriors(csmc.cliqSubFg)) > 0 || length(cliqSeparatorVarIds) > 1
    for si in cliqSeparatorVarIds
      vnd = getSolverData(getVariable(csmc.cliqSubFg, si), :parametric)
      beliefMsg.belief[si] = TreeBelief(deepcopy(vnd))
      # beliefMsg.belief[si] = TreeBelief(vnd.val, vnd.bw, vnd.inferdim, vnd.softtype.manifolds)
    end
  end

  for e in getEdgesParent(csmc.tree, csmc.cliq)
    @info "$(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)"
    putBeliefMessageUp!(csmc.tree, e, beliefMsg)# put!(csmc.tree.messages[e.index].upMsg, beliefMsg)
  end

  return waitForDown_ParametricStateMachine
end


"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 4
"""
function waitForDown_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-4, wait for down messages of needed")

  setCliqDrawColor(csmc.cliq, "turquoise")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  for e in getEdgesParent(csmc.tree, csmc.cliq)
    @info "$(csmc.cliq.index): take! on edge $(isa(e,Graphs.Edge) ? e.index : e)"
    # Blocks until data is available.
    beliefMsg = takeBeliefMessageDown!(csmc.tree, e) # take!(csmc.tree.messages[e.index].downMsg)
    @info "$(csmc.cliq.index): Belief message recieved with status $(beliefMsg.status)"


    #save down messages in msgsDown
    if beliefMsg.status == DOWNSOLVED
      setDwnMsg!(csmc, beliefMsg)
      # csmc.msgsDown = beliefMsg
      # push!(csmc.msgsDown, beliefMsg)

    else
      setCliqDrawColor(csmc.cliq, "red")
      # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

      @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
        @info "Par-4, $(csmc.cliq.index): put! error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        @async putBeliefMessageDown!(csmc.tree, e, LikelihoodMessage(ERROR_STATUS))#put!(csmc.tree.messages[e.index].downMsg,  LikelihoodMessage(ERROR_STATUS))
      end
      @error "Par-4, $(csmc.cliq.index): Exit with error state"
      return IncrementalInference.exitStateMachine
    end
  end

  return solveDown_ParametricStateMachine
end


"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 5
"""
function solveDown_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-5, Solving down")

  setCliqDrawColor(csmc.cliq, "red")
  # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  #TODO maybe change to symbols
  # for downmsgs in csmc.msgsDown
  downmsgs = getDwnMsgs(csmc) # csmc.msgsDown
    # TODO
    # updateMsgSeparators!(csmc.cliqSubFg, downmsgs)
    svars = getCliqSeparatorVarIds(csmc.cliq)
    # svars = DFG.listVariables(csmc.cliqSubFg)
    for (msym, belief) in downmsgs.belief
      if msym in svars
        #TODO maybe combine variable and factor in new prior?
        vnd = getSolverData(getVariable(csmc.cliqSubFg, msym), :parametric)
        @info "$(csmc.cliq.index): Updating separator $msym from message $(belief.val)"
        vnd.val .= belief.val
        vnd.bw .= belief.bw
      end
    end
  # end

  # store the cliqSubFg for later debugging
  # NOTE ITS not changed for now but keep here for possible future use
  # opts = getSolverParams(csmc.dfg)
  # if opts.dbg
  #   DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedownsolve"))
  #   drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforedownsolve.pdf"))
  # end

  #TODO DownSolve cliqSubFg
  #only down solve if its not a root
  if length(getParent(csmc.tree, csmc.cliq)) != 0
    frontals = getCliqFrontalVarIds(csmc.cliq)
    vardict, result, flatvars, Σ = solveConditionalsParametric(csmc.cliqSubFg, frontals)
    #TEMP testing difference
    # vardict, result = solveFactorGraphParametric(csmc.cliqSubFg)
    # Pack all results in variables
    if result.g_converged || result.f_converged
      @info "$(csmc.cliq.index): subfg optim converged updating variables"
      for (v,val) in vardict
        @info "$(csmc.cliq.index) down: updating $v : $val"
        vnd = getSolverData(getVariable(csmc.cliqSubFg, v), :parametric)
        #Update subfg variables
        vnd.val .= val.val
        vnd.bw .= val.cov
      end
    else
      @error "Par-5, clique $(csmc.cliq.index) failed to converge in down solve" result

      #propagate error to cleanly exit all cliques?
      beliefMsg = LikelihoodMessage(ERROR_STATUS)
      @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
        @info "Par-5, $(csmc.cliq.index): put! error on edge $(isa(e,Graphs.Edge) ? e.index : e)"
        @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messages[e.index].downMsg, beliefMsg)
      end
      @error "Par-5, $(csmc.cliq.index): Exit with error state"
      return IncrementalInference.exitStateMachine

    end
  end

  #TODO fill in belief
  cliqFrontalVarIds = getCliqFrontalVarIds(csmc.cliq)
  #TODO createBeliefMessageParametric
  # beliefMsg = createBeliefMessageParametric(csmc.cliqSubFg, cliqFrontalVarIds, solvekey=opts.solvekey)
  beliefMsg = LikelihoodMessage(DOWNSOLVED)
  for fi in cliqFrontalVarIds
    vnd = getSolverData(getVariable(csmc.cliqSubFg, fi), :parametric)
    beliefMsg.belief[fi] = TreeBelief(vnd)
    # beliefMsg.belief[fi] = TreeBelief(vnd.val, vnd.bw, vnd.inferdim, vnd.softtype.manifolds)
    @info "$(csmc.cliq.index): down message $fi : $beliefMsg"
  end

  #TODO sendBeliefMessageParametric(csmc, beliefMsg)
  #TODO maybe send a specific message to only the child that needs it
  @sync for e in getEdgesChildren(csmc.tree, csmc.cliq)
    @info "$(csmc.cliq.index): put! on edge $(isa(e,Graphs.Edge) ? e.index : e)"
    @async putBeliefMessageDown!(csmc.tree, e, beliefMsg)#put!(csmc.tree.messages[e.index].downMsg, beliefMsg)
  end

  @info "$(csmc.cliq.index): Solve completed"

  if isa(csmc.dfg, DFG.InMemoryDFGTypes)
    #TODO update frontal variables here directly
    frontsyms = getFrontals(csmc.cliq)
    transferUpdateSubGraphParametric!(csmc.dfg, csmc.cliqSubFg, frontsyms)
    #solve finished change color
    setCliqDrawColor(csmc.cliq, "lightblue")
    # csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

    @info "$(csmc.cliq.index): Finished"
    return IncrementalInference.exitStateMachine
  else
    #seems like a nice place to update remote variables here
    return updateRemote_ParametricStateMachine
  end
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 6
"""
function updateRemote_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-6, Updating Remote")
  #TODO update frontal variables remotely here
  #NOTE a new state is made to allow for coms retries and error traping kind of behaviour

  @info "$(csmc.cliq.index): Finish en klaar"
  return IncrementalInference.exitStateMachine

end


"""
    $SIGNATURES

Transfer contents of `src::AbstractDFG` variables `syms::Vector{Symbol}` to `dest::AbstractDFG`.

Notes
- Reads, `dest` := `src`, for all `syms`
"""
function transferUpdateSubGraphParametric!(dest::InMemoryDFGTypes,
                                           src::InMemoryDFGTypes,
                                           syms::Vector{Symbol},
                                           solveKey::Symbol=:parametric,
                                           logger=ConsoleLogger()  )
  #
  with_logger(logger) do
    @info "transferUpdateSubGraph! -- syms=$syms"
  end
  for v in getVariables(src)
    println("\n ", v.label,": ",  getSolverData(v, :parametric).val[1])
  end

  # #TEMP force the solver data
  # for v in syms
  #   getSolverData(getVariable(dest,v),:parametric).val .= getSolverData(getVariable(src, v),:parametric).val
  # end
  DFG.updateVariableSolverData!(dest, [getVariable(src, vIdx) for vIdx in syms], :parametric)
  #TODO this does not work
  # DFG.mergeUpdateGraphSolverData!(dest, src, syms)
  nothing
end
