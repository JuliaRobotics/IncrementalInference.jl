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

  # TODO JT toets om priors van seperators af te haal
  removedIds = removeSeperatorPriorsFromSubgraph!(csmc.cliqSubFg, csmc.cliq)
  infocsm(csmc, "Par-1 Removed ids $removedIds")


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

#TODO get this to work
struct ParametricMsgPrior{T} <: IncrementalInference.FunctorSingleton
  Z::Vector{Float64}
  inferdim::Float64
end

#TODO pack and unpack
getSample(s::ParametricMsgPrior, N::Int=1) = (s.Z, )

struct PackedParametricMsgPrior <: PackedInferenceType where T
  Z::String
  inferdim::Float64
end

function convert(::Type{PackedParametricMsgPrior}, d::ParametricMsgPrior)
  PackedParametricMsgPrior(string(d.Z), d.inferdim)
end
function convert(::Type{ParametricMsgPrior}, d::PackedParametricMsgPrior)
  ParametricMsgPrior([0.0], d.inferdim)
end

#TODO move to TreeBasedInitialization.jl
function addMsgFactors!(subfg::G,
                        msgs::ParametricBeliefMessage)::Vector{DFGFactor} where G <: AbstractDFG
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  svars = DFG.getVariableIds(subfg)
  for (msym, belief) = (msgs.belief)
    if msym in svars
      # TODO prior missing manifold information? is it not available in variable?
      # TODO add MsgPrior parametric type
      @warn "TODO I just packed in MvNormal and ignored cov"
      msgPrior =  MsgPrior(MvNormal(belief.val[:,1], 1), belief.inferdim)
      fc = addFactor!(subfg, [msym], msgPrior, autoinit=false)
      push!(msgfcts, fc)
    end
  end
  return msgfcts
end

"""
    $SIGNATURES

Notes
- Parametric state machine function nr. 2
"""
function waitForUp_ParametricStateMachine(csmc::CliqStateMachineContainer)

  infocsm(csmc, "Par-2, wait for up messages of needed")

  setCliqDrawColor(csmc.cliq, "purple") #TODO don't know if this is correct color
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  for e in Graphs.out_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): take! on edge $(e.index)"
    # Blocks until data is available.
    beliefMsg = take!(csmc.tree.messages[e.index].upMsg)
    @info "$(csmc.cliq.index): Belief message recieved with status $(beliefMsg.status)"
    #TODO save up message (and add priors to cliqSubFg) gaan eers in volgende stap kyk hoe dit werk
    #kies csmc vir boodskappe vir debugging, dis 'n vector een per kind knoop
    push!(csmc.parametricMsgsUp, beliefMsg)
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

  setCliqDrawColor(csmc.cliq, "red")
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  # is dit nodig om vir parametric die factors by te sit?
  #TODO maybe change to symbols
  msgfcts = DFGFactor[]
  for upmsgs in csmc.parametricMsgsUp
    append!(msgfcts, addMsgFactors!(csmc.cliqSubFg, upmsgs))
  end

  # store the cliqSubFg for later debugging
  opts = getSolverParams(csmc.dfg)
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_beforeupsolve.pdf"))
  end

  #TODO UpSolve cliqSubFg here
  # slaap om somme te simileer
  sleep(rand()*2)
  # solveFactorGraphParametric!(csmc.cliqSubFg)

  # Done with solve delete factors
  deleteMsgFactors!(csmc.cliqSubFg, msgfcts)

  # store the cliqSubFg for later debugging
  if opts.dbg
    DFG.saveDFG(csmc.cliqSubFg, joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve"))
    drawGraph(csmc.cliqSubFg, show=false, filepath=joinpath(opts.logpath,"logs/cliq$(csmc.cliq.index)/fg_afterupsolve.pdf"))
  end

  #TODO fill in belief
  # createBeliefMessageParametric(csmc.cliqSubFg, csmc.cliq, solvekey=opts.solvekey)
  cliqSeparatorVarIds = getCliqSeparatorVarIds(csmc.cliq)
  beliefMsg = ParametricBeliefMessage(upsolved)
  for si in cliqSeparatorVarIds
    vnd = getVariableData(csmc.cliqSubFg, si, solveKey=:parametric)
    beliefMsg.belief[si] = TreeBelief(vnd.val, vnd.bw, vnd.inferdim)
  end

  #NOTE Graphs.jl in_edges does not work. So extended above
  for e in in_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): put! on edge $(e.index)"
    put!(csmc.tree.messages[e.index].upMsg, beliefMsg)
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

  setCliqDrawColor(csmc.cliq, "turquoise")
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  for e in Graphs.in_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): take! on edge $(e.index)"
    # Blocks until data is available.
    beliefMsg = take!(csmc.tree.messages[e.index].downMsg)
    @info "$(csmc.cliq.index): Belief message recieved with status $(beliefMsg.status)"
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

  setCliqDrawColor(csmc.cliq, "red")
  csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

  #TODO UpDown cliqSubFg here
  # slaap om somme te simileer
  sleep(rand()*2)

  #TODO fill in belief
  beliefMsg = ParametricBeliefMessage(downsolved)

  #TODO send a specific message to only the child that needs it
  for e in out_edges(csmc.cliq, csmc.tree.bt)
    @info "$(csmc.cliq.index): put! on edge $(e.index)"
    put!(csmc.tree.messages[e.index].downMsg, beliefMsg)
  end

  @info "$(csmc.cliq.index): Solve Finished"

  if isa(csmc.dfg, DFG.InMemoryDFGTypes)
    #TODO update frontal variables here directly

    #solve finished change color
    setCliqDrawColor(csmc.cliq, "lightblue")
    csmc.drawtree ? drawTree(csmc.tree, show=false, filepath=joinpath(getSolverParams(csmc.dfg).logpath,"bt.pdf")) : nothing

    @info "$(csmc.cliq.index): Finish en klaar"
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

Transfer contents of `src::FactorGraph` variables `syms::Vector{Symbol}` to `dest::FactorGraph`.

Notes
- Reads, `dest` := `src`, for all `syms`
"""
function transferUpdateSubGraphParametric!(dest::T,
                                           src::T,
                                           syms::Vector{Symbol},
                                           solveKey::Symbol=:parametric,
                                           logger=ConsoleLogger()  ) where {T <: DFG.InMemoryDFGTypes}
  #
  with_logger(logger) do
    @info "transferUpdateSubGraph! -- syms=$syms"

    DFG.mergeUpdateVariableSolverData(dest, src, syms)

  end
  nothing
end
