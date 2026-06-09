
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
  bt::MetaDiGraph{Int, Float64}
  btid::Int
  frontals::Dict{Symbol, CliqueId{Int}}
  eliminationOrder::Vector{Symbol}
  buildTime::Float64
end

const BayesTree = MetaBayesTree


@kwdef mutable struct CSMOptions
  solverparams::Any #FIXME deprecation step
  verbose::Bool = false
  verbosefid::Any = stdout
  drawtree::Bool = false
  limititers::Int = -1
  downsolve::Bool = false
  upsolve::Bool = true # TODO unchecked consolidation, consolidate from solverparams
  incremental::Bool = false
  solve_progressbar::Any = nothing
  algorithm::Symbol = :default
  solveKey::Symbol = algorithm
  recordhistory::Bool = false
  delay::Bool = false
  multithread::Bool = false
  timeout::Union{Nothing, <:Real} = nothing
  delaycliqs::Vector{Symbol} = Symbol[]
  recordcliqs::Vector{Symbol} = Symbol[]
  skipcliqids::Vector{Symbol} = Symbol[]
  limititercliqs::Vector{Pair{Symbol, Int}} = Pair{Symbol, Int}[]
end


"""
    $TYPEDEF

Container for upward tree solve / initialization.

DevNotes
- TODO more direct clique access (cliq, parent, children), for multi-process solves
"""
@kwdef mutable struct CliqStateMachineContainer{
  BTND,
  G <: AbstractDFG,
  InMemG <: InMemoryDFGTypes,
  BT <: AbstractBayesTree,
}
  dfg::G
  cliqSubFg::InMemG
  tree::BT
  cliq::TreeClique
    incremental::Bool
    drawtree::Bool
  dodownsolve::Bool
    delay::Bool
  opts::SolverParams
  refactoring::Dict{Symbol, String} = Dict{Symbol, String}()
  oldcliqdata::BTND = BayesTreeNodeData()
  logger::SimpleLogger = SimpleLogger(Base.stdout)
  cliqId::CliqueId = cliq.id # obsolete?
    algorithm::Symbol = :default
  init_iter::Int = 0
  enableLogging::Bool = true
    solveKey::Symbol = :default
  _csm_iter::Int = 0
end


#TODO use @NamedTuple if julia compat > 1.5

const CSMHistoryTuple = NamedTuple{
  (:timestamp, :id, :f, :csmc),
  Tuple{DateTime, Int, Function, CliqStateMachineContainer},
}
const CSMHistory = Vector{CSMHistoryTuple}

#
