# STATUS messages [:initialized;:upsolved;:marginalized;:downsolved;:uprecycled]

"""
    $SIGNATURES

EXPERIMENTAL: Init and start state machine for parametric solve.

Notes:
- will call on values from children or parent cliques
- can be called multiple times
- Assumes all cliques in tree are being solved simultaneously and in similar manner.
- State machine rev.1 -- copied from first TreeBasedInitialization.jl.
- Doesn't do partial initialized state properly yet.
"""
function initStartCliqStateMachineParametric!(dfg::G,
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

  #csmc = CliqStateMachineContainer(dfg, initfg(destType), tree, cliq, prnt, children, false, incremental, drawtree, downsolve, delay, getSolverParams(dfg), oldcliqdata, logger)
  csmc = CliqStateMachineContainer(dfg, initfg(destType, params=getSolverParams(dfg)), tree, cliq, prnt, children, false, incremental, drawtree, downsolve, delay, getSolverParams(dfg), oldcliqdata, logger)

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
  # NOTE add all frontal factor neighbors DEV CASE -- use getData(cliq).dwnPotentials instead
  # fnsyms = getCliqVarsWithFrontalNeighbors(csmc.dfg, csmc.cliq)

  infocsm(csmc, "Par-1, build subgraph syms=$(syms)")

  buildSubgraphFromLabels!(csmc.dfg, csmc.cliqSubFg, syms)
  # TODO review, are all updates atomic???
  # if isa(csmc.dfg, DFG.InMemoryDFGTypes)
  #   csmc.cliqSubFg = csmc.dfg
  # else
  #  buildSubgraphFromLabels!(dfg, csmc.cliqSubFg, syms)
  # end

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

  for e in Graphs.out_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): take! on edge $(e.index)"
    # Blocks until data is available.
    believeMsg = take!(csmc.tree.messages[e.index].upMsg)
    @info "$(csmc.cliq.index): Believe message recieved with status $(believeMsg.status)"
    #TODO save up message and add priors to cliqSubFg
  end

  return solveUp_ParametricStateMachine
end

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

  #TODO UpSolve cliqSubFg here
  # slaap om somme te simileer
  sleep(rand()*2)

  #TODO fill in belief
  believeMsg = ParametricBelieveMessage(:upsolved)

  #NOTE Graphs.jl in_edges does not work. So extended above
  for e in in_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): put! on edge $(e.index)"
    put!(csmc.tree.messages[e.index].upMsg, believeMsg)
  end

  return waitForDown_ParametricStateMachine
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 4
"""
function waitForDown_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-4, wait for up messages of needed")

  for e in Graphs.in_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): take! on edge $(e.index)"
    # Blocks until data is available.
    believeMsg = take!(csmc.tree.messages[e.index].downMsg)
    @info "$(csmc.cliq.index): Believe message recieved with status $(believeMsg.status)"
    #TODO save up message and add priors to cliqSubFg
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


  #TODO UpDown cliqSubFg here
  # slaap om somme te simileer
  sleep(rand()*2)

  #TODO fill in belief
  believeMsg = ParametricBelieveMessage(:downSolved)

  #TODO send a specific message to only the child that needs it
  for e in out_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): put! on edge $(e.index)"
    put!(csmc.tree.messages[e.index].downMsg, believeMsg)
  end


  if isa(csmc.dfg, DFG.InMemoryDFGTypes)
    #in memory type exit as variables should be up to date
    #Finish en klaar
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

  #Finish en klaar
  return IncrementalInference.exitStateMachine

end
