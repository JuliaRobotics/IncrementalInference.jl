
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
  solve_progressbar::Any = nothing
  # execution options
  solverparams::Any #FIXME deprecation step
  downsolve::Bool = false
  upsolve::Bool = true # TODO unchecked consolidation, consolidate from solverparams
  incremental::Bool = false
  algorithm::Symbol = :default
  solveKey::Symbol = algorithm
  delay::Bool = false
  multithread::Bool = false
  # debug options
  limititers::Int = -1
  recordhistory::Bool = false
  storeOld::Bool = false
  verbose::Bool = false
  verbosefid::Any = stdout
  drawtree::Bool = false
  dotreedraw::Vector{Int} = Int[1;]
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
  dodownsolve::Bool = false
  opts::SolverParams = getSolverParams(cliqSubFg)
  refactoring::Dict{Symbol, String} = Dict{Symbol, String}()
  oldcliqdata::BTND = BayesTreeNodeData()
  logger::SimpleLogger = SimpleLogger(Base.stdout)
  cliqId::CliqueId = cliq.id # obsolete?
  init_iter::Int = 0
  enableLogging::Bool = true
  _csm_iter::Int = 0
  csmoptions::CSMOptions = CSMOptions(solverparams = opts)
end


#TODO use @NamedTuple if julia compat > 1.5

const CSMHistoryTuple = NamedTuple{
  (:timestamp, :id, :f, :csmc),
  Tuple{DateTime, Int, Function, CliqStateMachineContainer},
}
const CSMHistory = Vector{CSMHistoryTuple}

#
