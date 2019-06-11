

# BayesTree declarations
"""
$(TYPEDEF)

Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::FactorGraph`.
"""
mutable struct BayesTree
  bt
  btid::Int
  cliques::Dict{Int,Graphs.ExVertex}
  frontals::Dict{Symbol,Int}
end

function emptyBayesTree()
    bt =   BayesTree(Graphs.inclist(Graphs.ExVertex,is_directed=true),
                     0,
                     Dict{Int,Graphs.ExVertex}(),
                     #[],
                     Dict{AbstractString, Int}())
    return bt
end

"""
    $TYPEDEF

Container for upward tree solve / initialization.

TODO
- remove proceed
- more direct clique access (cliq, parent, children), for multi-process solves
"""
mutable struct CliqStateMachineContainer
  dfg::DistributedFactorGraphs.GraphsDFG
  cliqSubFg::DistributedFactorGraphs.GraphsDFG
  tree::BayesTree
  cliq::Graphs.ExVertex
  parentCliq::Vector{Graphs.ExVertex}
  childCliqs::Vector{Graphs.ExVertex}
  forceproceed::Bool # TODO: bad flag that must be removed by refactoring sm
  incremental::Bool
  drawtree::Bool
  dodownsolve::Bool
  refactoring::Dict{Symbol, String}
  CliqStateMachineContainer() = new()
  CliqStateMachineContainer(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) = new(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,Dict{Symbol,String}())
end

"""
$(TYPEDEF)

Data structure for each clique in the Bayes (Junction) tree.
"""
mutable struct BayesTreeNodeData
  frontalIDs::Vector{Symbol}
  conditIDs::Vector{Symbol}
  inmsgIDs::Vector{Symbol} # Int
  potIDs::Vector{Symbol} # Int # this is likely redundant TODO -- remove
  potentials::Vector{Symbol}
  partialpotential::Vector{Bool}
  cliqAssocMat::Array{Bool,2}
  cliqMsgMat::Array{Bool,2}
  directvarIDs::Vector{Symbol} # Int
  directFrtlMsgIDs::Vector{Symbol} # Int
  msgskipIDs::Vector{Symbol} # Int
  itervarIDs::Vector{Symbol} # Int
  directPriorMsgIDs::Vector{Symbol} # Int
  debug
  debugDwn
  # future might concentrate these four fields down to two
  upMsg::Dict{Symbol, BallTreeDensity}
  dwnMsg::Dict{Symbol, BallTreeDensity}
  upInitMsgs::Dict{Int, Dict{Symbol, BallTreeDensity}}
  downInitMsg::Dict{Symbol, BallTreeDensity}

  allmarginalized::Bool
  initialized::Symbol
  upsolved::Bool
  downsolved::Bool
  initUpChannel::Channel{Symbol}
  initDownChannel::Channel{Symbol}
  solveCondition::Condition
  statehistory::Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}
  BayesTreeNodeData() = new()
  BayesTreeNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],
                                x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],
                                x[21], x[22], x[23], x[24], x[25], x[26],
                                Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}() )
end

# TODO -- this should be a constructor
function emptyBTNodeData()
  BayesTreeNodeData(Int[],Int[],Int[],
                    Int[],Int[],Bool[],
                    Array{Bool}(undef, 0,0),
                    Array{Bool}(undef, 0,0),
                    Int[],Int[],
                    Int[],Int[],Int[],
                    nothing, nothing,
                    Dict{Symbol, BallTreeDensity}(),  # :null => AMP.manikde!(zeros(1,1), [1.0;], (:Euclid,))),
                    Dict{Symbol, BallTreeDensity}(),  # :null => AMP.manikde!(zeros(1,1), [1.0;], (:Euclid,))),
                    Dict{Int, Dict{Symbol, BallTreeDensity}}(),
                    Dict{Symbol, BallTreeDensity}(),
                    false, :null,
                    false, false,
                    Channel{Symbol}(1), Channel{Symbol}(1), Condition()  )
end



"""
$(TYPEDEF)
"""
mutable struct NBPMessage <: Singleton
  p::Dict{Symbol, EasyMessage}
end

"""
$(TYPEDEF)
"""
mutable struct PotProd
    Xi::Symbol # Int
    prev::Array{Float64,2}
    product::Array{Float64,2}
    potentials::Array{BallTreeDensity,1}
    potentialfac::Vector{Symbol}
end
"""
$(TYPEDEF)
"""
mutable struct CliqGibbsMC
    prods::Array{PotProd,1}
    lbls::Vector{Symbol}
    CliqGibbsMC() = new()
    CliqGibbsMC(a,b) = new(a,b)
end
"""
$(TYPEDEF)
"""
mutable struct DebugCliqMCMC
    mcmc::Union{Nothing, Array{CliqGibbsMC,1}}
    outmsg::NBPMessage
    outmsglbls::Dict{Symbol, Symbol} # Int
    priorprods::Vector{CliqGibbsMC}
    DebugCliqMCMC() = new()
    DebugCliqMCMC(a,b,c,d) = new(a,b,c,d)
end

"""
$(TYPEDEF)
"""
mutable struct UpReturnBPType
  upMsgs::NBPMessage
  dbgUp::DebugCliqMCMC
  IDvals::Dict{Symbol, EasyMessage} # Int
  keepupmsgs::Dict{Symbol, BallTreeDensity} # TODO Why separate upMsgs?
  totalsolve::Bool
  UpReturnBPType() = new()
  UpReturnBPType(x1,x2,x3,x4,x5) = new(x1,x2,x3,x4,x5)
end

"""
$(TYPEDEF)
"""
mutable struct DownReturnBPType
    dwnMsg::NBPMessage
    dbgDwn::DebugCliqMCMC
    IDvals::Dict{Symbol,EasyMessage} # Int
    keepdwnmsgs::Dict{Symbol, BallTreeDensity}
end

"""
$(TYPEDEF)
"""
mutable struct FullExploreTreeType{T, T2}
  fg::DistributedFactorGraphs.GraphsDFG
  bt::T2
  cliq::Graphs.ExVertex
  prnt::T
  sendmsgs::Vector{NBPMessage}
end

const ExploreTreeType{T} = FullExploreTreeType{T, BayesTree}
const ExploreTreeTypeLight{T} = FullExploreTreeType{T, Nothing}


function ExploreTreeType(fgl::G,
                         btl::BayesTree,
                         vertl::Graphs.ExVertex,
                         prt::T,
                         msgs::Array{NBPMessage,1} ) where {G <: AbstractDFG, T}
  #
  ExploreTreeType{T}(fgl, btl, vertl, prt, msgs)
end

"""
$(TYPEDEF)
"""
mutable struct MsgPassType
  fg::GraphsDFG
  cliq::Graphs.ExVertex
  vid::Symbol # Int
  msgs::Array{NBPMessage,1}
  N::Int
end
