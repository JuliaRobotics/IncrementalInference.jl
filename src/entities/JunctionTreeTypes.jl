

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

# NOTE Temporary fix? Overwrite deepcopy on MetaBayesTree to strip out copying the channels. 
# see https://github.com/JuliaRobotics/IncrementalInference.jl/issues/1530
# possible fix in https://github.com/JuliaLang/julia/pull/46406
import Base.deepcopy_internal
function Base.deepcopy_internal(bt::MetaBayesTree, stackdict::IdDict)
 
  if haskey(stackdict, bt)
      return stackdict[bt]
  end

  mg = bt.bt
  
  graph = deepcopy_internal(mg.graph, stackdict)
  vprops = deepcopy_internal(mg.vprops, stackdict)
  T = eltype(mg)
  # dropping all edge data
  eprops = Dict{MetaGraphs.SimpleEdge{T},MetaGraphs.PropDict}()
  gprops = deepcopy_internal(mg.gprops, stackdict)
  weightfield = deepcopy_internal(mg.weightfield, stackdict)
  defaultweight = deepcopy_internal(mg.defaultweight, stackdict)
  metaindex = deepcopy_internal(mg.metaindex, stackdict)
  indices = deepcopy_internal(mg.indices, stackdict)

  mg_cpy = MetaDiGraph(graph, vprops, eprops, gprops, weightfield, defaultweight, metaindex, indices)

  bt_cpy = MetaBayesTree(mg_cpy, bt.btid, deepcopy_internal(bt.frontals, stackdict), deepcopy_internal(bt.eliminationOrder, stackdict), bt.buildTime)

  stackdict[bt] = bt_cpy
  return bt_cpy
end

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


#