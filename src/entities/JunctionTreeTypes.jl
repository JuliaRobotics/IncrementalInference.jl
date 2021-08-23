
export MetaBayesTree, BayesTree
export CSMHistoryTuple

## ========================================================================================================================
## Bayes Trees
## ========================================================================================================================

abstract type AbstractBayesTree end


# TODO DEV MetaGraphs bayes tree, will potentially also make a LightBayesTree, CloudBayesTree,
"""
$(TYPEDEF)
Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::AbstractDFG`.
"""
mutable struct MetaBayesTree <: AbstractBayesTree
  bt::MetaDiGraph{Int,Float64}
  btid::Int
  frontals::Dict{Symbol,CliqueId{Int}}
  eliminationOrder::Vector{Symbol}
  buildTime::Float64
end

const BayesTree = MetaBayesTree

MetaBayesTree() = MetaBayesTree(MetaDiGraph{Int,Float64}(), 0, Dict{AbstractString, Int}(), Symbol[], 0.0)

Base.propertynames(x::MetaBayesTree, private::Bool=false) = (:bt, :btid, :cliques, :frontals, :variableOrder, :buildTime)

Base.getproperty(x::MetaBayesTree,f::Symbol) = begin
    if f == :cliques
      @warn "Don't use cliques field directly, use eg. getClique(tree, cliqId)" maxlog=1
      d = Dict{Int,Any}()
      for (k,v) in x.bt.vprops
        d[k] = v[:clique]
      end
      return d
    elseif f == :variableOrder
      @warn "MetaBayesTree.variableOrder has been deprecated in favor of .eliminationOrder"
      return getfield(x,:eliminationOrder)
    else
      getfield(x,f)
    end
  end

function Base.setproperty!(x::MetaBayesTree, f::Symbol, val)
  if f == :cliques
    @warn "`setproperty!(clique)` Don't use cliques field directly, use eg. addClique(tree, cliqId)" maxlog=1
    for (k,v) in val
      set_prop!(x.bt, k, :clique, v)
    end
  else
    setfield!(x,f,val)
  end
end

function getMessageChannels(tree::MetaBayesTree)

  d = Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{LikelihoodMessage},Channel{LikelihoodMessage}}}}()

  for (k,e) in tree.bt.eprops
    d[k.dst] = (upMsg=e[:upMsg], downMsg=e[:downMsg])
  end
  
  return d

end

function Base.show(io::IO, mbt::MetaBayesTree)
  printstyled(io, "MetaBayesTree\n", color=:blue)
  println(io, "  Nr cliques:  ", length(mbt.cliques))

  # TODO ad dmore stats: max depth, widest point, longest chain, max clique size, average nr children

  nothing
end

Base.show(io::IO, ::MIME"text/plain", mbt::MetaBayesTree) = show(io, mbt)


"""
    $TYPEDEF

Container for upward tree solve / initialization.

DevNotes
- TODO more direct clique access (cliq, parent, children), for multi-process solves
"""
mutable struct CliqStateMachineContainer{BTND, G <: AbstractDFG, InMemG <: InMemoryDFGTypes, BT <: AbstractBayesTree}
  dfg::G
  cliqSubFg::InMemG
  tree::BT
  cliq::TreeClique
  incremental::Bool
  drawtree::Bool
  dodownsolve::Bool
  delay::Bool
  opts::SolverParams
  refactoring::Dict{Symbol, String}
  oldcliqdata::BTND
  logger::SimpleLogger
  cliqId::CliqueId
  algorithm::Symbol
  init_iter::Int
  enableLogging::Bool
  solveKey::Symbol
  _csm_iter::Int
end

#TODO use @NamedTuple if julia compat > 1.5

const CSMHistoryTuple =  NamedTuple{(:timestamp, :id, :f, :csmc), Tuple{DateTime, Int, Function, CliqStateMachineContainer}}
const CSMHistory = Vector{CSMHistoryTuple}

Base.show(io::IO, o::CSMHistoryTuple) = print(io, "$(o[1]), $(o[2]), $(o[3])")


function CliqStateMachineContainer( dfg::G,
                                    cliqSubFg::M,
                                    tree::T,
                                    cliq::TreeClique,
                                    incremental::Bool,
                                    drawtree::Bool,
                                    dodownsolve::Bool,
                                    delay::Bool,
                                    opts::SolverParams,
                                    refactoring::Dict{Symbol,String}=Dict{Symbol,String}(),
                                    oldcliqdata::BTND=BayesTreeNodeData(),
                                    logger::SimpleLogger=SimpleLogger(Base.stdout);
                                    cliqId::CliqueId = cliq.id,
                                    algorithm::Symbol = :default,
                                    init_iter::Int=0,
                                    enableLogging::Bool=true,
                                    solveKey::Symbol=:default,
                                    _csm_iter::Int=0  ) where {BTND, G <: AbstractDFG, M <: InMemoryDFGTypes, T <: AbstractBayesTree}
  #
  CliqStateMachineContainer{BTND, G, M, T}( dfg,
                                            cliqSubFg,
                                            tree,
                                            cliq,
                                            incremental,
                                            drawtree,
                                            dodownsolve,
                                            delay,
                                            opts,
                                            refactoring,
                                            oldcliqdata,
                                            logger,
                                            cliqId,
                                            algorithm,
                                            init_iter,
                                            enableLogging,
                                            solveKey,
                                            _csm_iter )
  #
end

# TODO resolve name conflict
function DFG.compare( cs1::CliqStateMachineContainer{BTND1, T1, InMemG1, BT1},
                      cs2::CliqStateMachineContainer{BTND2, T2, InMemG2, BT2};
                      skip::Vector{Symbol}=Symbol[] ) where {BTND1, T1 <: AbstractDFG, InMemG1 <: InMemoryDFGTypes, BT1 <: AbstractBayesTree, BTND2, T2 <: AbstractDFG, InMemG2 <: InMemoryDFGTypes, BT2 <: AbstractBayesTree}
  #
  BTND1 == BTND2 ? nothing : @warn("oldcliqdata::$BTND1 != oldcliqdata::$BTND2")
  T1 == T2 ? nothing : @warn("dfg::$T1 != dfg::$T2")
  InMemG1 == InMemG2 ? nothing : @warn("cliqSubFg::$InMemG1 != cliqSubFg::$InMemG2")
  BT1 == BT2 ? nothing : @warn("tree::$BQ1 != tree::$BT2")

  TP = true
  @warn "Skipping compare of CSMC.dfg and .cliqSubFg"
  # TP = TP && compare(cs1.dfg,  cs2.dfg)
  # TP = TP && compare(cs1.cliqSubFg,  cs2.cliqSubFg)
  @warn "Skipping compare of CSMC.tree"
  # TP = TP && compare(cs1.tree,  cs2.tree)
  TP = TP && compare(cs1.cliq,  cs2.cliq)
  TP = TP && compare(cs1.cliqId,  cs2.cliqId)
  TP = TP && length(cs1.parentCliq) == length(cs2.parentCliq)
  for i in 1:length(cs1.parentCliq)
    TP = TP && compare(cs1.parentCliq[i],  cs2.parentCliq[i])
  end
  TP = TP && length(cs1.childCliqs) == length(cs2.childCliqs)
  for i in 1:length(cs1.childCliqs)
    TP = TP && compare(cs1.childCliqs[i],  cs2.childCliqs[i])
  end
  TP = TP && compare(cs1.incremental,  cs2.incremental)
  TP = TP && compare(cs1.drawtree,  cs2.drawtree)
  TP = TP && compare(cs1.dodownsolve,  cs2.dodownsolve)
  TP = TP && compare(cs1.delay,  cs2.delay)
  @warn "skipping compare on csmc.opts::SolverParams"
  # TP = TP && compare(cs1.opts,  cs2.opts)
  TP = TP && compare(cs1.refactoring,  cs2.refactoring)
  # TP = TP && compare(cs1.oldcliqdata,  cs2.oldcliqdata)
  # TP = TP && compare(cs1.logger,  cs2.logger)

  return TP
end


#