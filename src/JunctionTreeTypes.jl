

## Bayes Trees

abstract type AbstractBayesTree end

#TODO - see #540 related to indexing and ids

# BayesTree declarations
const BTGdict = GenericIncidenceList{TreeClique,Edge{TreeClique},Array{TreeClique,1},Array{Array{Edge{TreeClique},1},1}}

"""
$(TYPEDEF)

Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::AbstractDFG`.
"""
mutable struct BayesTree <: AbstractBayesTree
  bt::BTGdict
  btid::Int
  cliques::Dict{Int,TreeClique}
  frontals::Dict{Symbol,Int}
  #TEMP JT for evaluation, store message channels associated with edges between nodes Int -> edge id.  TODO rather store in graph
  messages::Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{LikelihoodMessage},Channel{LikelihoodMessage}}}}
  variableOrder::Vector{Symbol}
  buildTime::Float64
end

BayesTree() = BayesTree(Graphs.inclist(TreeClique,is_directed=true),
                        0,
                        Dict{Int,TreeClique}(),
                        Dict{AbstractString, Int}(),
                        Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{LikelihoodMessage},Channel{LikelihoodMessage}}}}(),
                        Symbol[],
                        0.0  )

#NOTE select type for development
# emptyBayesTree() = BayesTree()
# emptyBayesTree() = MetaBayesTree()

#TEMP switch the tree to use NOTE under development don't use MetaBayesTree yet
global UseMetaBayesTree = false
setUseMetaBayesTree(b::Bool) = global UseMetaBayesTree = b
function emptyBayesTree()
  global UseMetaBayesTree
  if UseMetaBayesTree
    @warn "Experimental, do not use yet, MetaBayesTree is under development"
    return MetaBayesTree()
  else
    return BayesTree()
  end
end

# TODO DEV MetaGraphs bayes tree, will potentially also make a LightBayesTree, CloudBayesTree,
"""
$(TYPEDEF)
Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::AbstractDFG`.
"""
mutable struct MetaBayesTree <: AbstractBayesTree
  bt::MetaDiGraph{Int,Float64}
  btid::Int
  # cliques::Dict{Int,TreeClique}
  frontals::Dict{Symbol,Int}
  variableOrder::Vector{Symbol}
  buildTime::Float64
end

MetaBayesTree() = MetaBayesTree(MetaDiGraph{Int,Float64}(), 0, Dict{AbstractString, Int}(), Symbol[], 0.0)

Base.propertynames(x::MetaBayesTree, private::Bool=false) = (:bt, :btid, :cliques, :frontals, :variableOrder, :buildTime)

Base.getproperty(x::MetaBayesTree,f::Symbol) = begin
    if f == :cliques
      if !(@isdefined getCliquesWarnOnce)
        @warn "Maybe don't use cliques field directly, TODO implement add/update/get/delete eg. getClique(tree, cliqId)"
        global getCliquesWarnOnce = true
      end
      d = Dict{Int,Any}()
      for (k,v) in x.bt.vprops
        d[k] = v[:clique]
      end
      return d
    else
      getfield(x,f)
    end
  end

function Base.setproperty!(x::MetaBayesTree, f::Symbol, val)
  if f == :cliques
    if !(@isdefined setCliquesWarnOnce)
      @warn "Maybe don't use cliques field directly, TODO implement add/update/get/delete eg. getClique(tree, cliqId)"
      global setCliquesWarnOnce = true
    end
    for (k,v) in val
      set_prop!(x.bt, k, :clique, v)
    end
  else
    setfield!(x,f,val)
  end
end

function MetaBayesTree(tree::BayesTree)
  # create graph from Graphs.jl adjacency_matrix
  mtree = MetaBayesTree(MetaDiGraph{Int, Float64}(MetaGraphs.SimpleDiGraph(Graphs.adjacency_matrix(tree.bt))), tree.btid, tree.frontals, tree.variableOrder, tree.buildTime)

  #deep copy over properties
  for v in tree.bt.vertices
    # set_prop!(mtree.bt, v.index, :label, deepcopy(v.label))
    set_prop!(mtree.bt, v.index, :clique, deepcopy(v))
  end

  ##  FIXME: Use common location for channels #675 (DF, asking for BTND)
  ##  TODO: placeholder for edge stored Channels
  ## set message passing properties,
  # for e in MetaGraphs.edges(mtree.bt)
  #   set_prop!(mtree.bt, e, :upMsg, Channel{BelieveMessage}(0))
  #   set_prop!(mtree.bt, e, :downMsg, Channel{BelieveMessage}(0))
  # end

  return mtree

end

# A replacement for to_dot that saves only plotting attributes
function savedot_attributes(io::IO, g::MetaDiGraph)
    write(io, "digraph G {\n")
    for p in props(g)
        write(io, "$(p[1])=$(p[2]);\n")
    end

    for v in MetaGraphs.vertices(g)
        write(io, "$v")
        if length(props(g, v)) > 0
            write(io, " [ ")
        end
        for p in props(g, v)
            # key = p[1]
            # write(io, "$key=\"$(p[2])\",")
            for (k,v) in p[2]
              write(io, "\"$k\"=\"$v\",")
            end
        end
        if length(props(g, v)) > 0
            write(io, "];")
        end
        write(io, "\n")
    end

    for e in MetaGraphs.edges(g)
        write(io, "$(MetaGraphs.src(e)) -> $(MetaGraphs.dst(e)) [ ")
        if MetaGraphs.has_prop(g, e, :downMsg) && MetaGraphs.has_prop(g, e, :upMsg)
          if isready(MetaGraphs.get_prop(g, e, :downMsg))
            write(io, "color=red")
          elseif isready(MetaGraphs.get_prop(g, e, :upMsg))
            write(io, "color=orange")
          else
            write(io, "color=black")
          end
        end
        write(io, "]\n")
    end
    write(io, "}\n")
end

function Graphs.to_dot(mdigraph::MetaDiGraph)
  g = deepcopy(mdigraph)
  for (i,val) in g.vprops
    push!(g.vprops[i],:attributes=>val[:clique].attributes)
    delete!(g.vprops[i],:clique)
    delete!(g.vprops[i],:index)
  end
  m = PipeBuffer()
  savedot_attributes(m, g)
  data = take!(m)
  close(m)
  return String(data)
end

"""
    $TYPEDEF

Container for upward tree solve / initialization.

TODO
- remove proceed
- more direct clique access (cliq, parent, children), for multi-process solves
"""
mutable struct CliqStateMachineContainer{BTND, G <: AbstractDFG, InMemG <: InMemoryDFGTypes, BT <: AbstractBayesTree}
  dfg::G
  cliqSubFg::InMemG
  tree::BT
  cliq::TreeClique
  parentCliq::Vector{TreeClique}
  childCliqs::Vector{TreeClique}
  incremental::Bool
  drawtree::Bool
  dodownsolve::Bool
  delay::Bool
  opts::SolverParams
  refactoring::Dict{Symbol, String}
  oldcliqdata::BTND
  logger::SimpleLogger
  cliqKey::Int
  #TODO towards consolidated messages
  # Decision is pull/fetch-model #674 -- i.e. CSM only works inside its own csmc and fetchs messages from neighbors
  # NOTE, DF going for Dict over Vector
  # msgsUp::Dict{Int, LikelihoodMessage} # Vector{LikelihoodMessage}
  # msgsDown::LikelihoodMessage # Vector{LikelihoodMessage}
end

const CSMHistory = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}


function CliqStateMachineContainer( dfg::G,
                                    cliqSubFg::M,
                                    tree::T,
                                    cliq::TreeClique,
                                    parentCliq::Vector{TreeClique},
                                    childCliqs::Vector{TreeClique},
                                    incremental::Bool,
                                    drawtree::Bool,
                                    dodownsolve::Bool,
                                    delay::Bool,
                                    opts::SolverParams,
                                    refactoring::Dict{Symbol,String}=Dict{Symbol,String}(),
                                    oldcliqdata::BTND=BayesTreeNodeData(),
                                    logger::SimpleLogger=SimpleLogger(Base.stdout);
                                    cliqKey::Int = cliq.index) where {BTND, G <: AbstractDFG, M <: InMemoryDFGTypes, T <: AbstractBayesTree}
  #
  CliqStateMachineContainer{BTND, G, M, T}( dfg,
                                            cliqSubFg,
                                            tree,
                                            cliq,
                                            parentCliq,
                                            childCliqs,
                                            incremental,
                                            drawtree,
                                            dodownsolve,
                                            delay,
                                            opts,
                                            refactoring,
                                            oldcliqdata,
                                            logger,
                                            cliqKey )
  #
end


function compare(cs1::CliqStateMachineContainer{BTND1, T1, InMemG1, BT1},
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
  TP = TP && compare(cs1.cliqKey,  cs2.cliqKey)
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


## === DF THIS MUST BE DELETED
# JT Using this to fix parametric tree again 
# include additional Tx buffers for later use
# fech vs take! #855
# paremetric needs take! to work
mutable struct MessageStore
  upRx::Dict{Int, LikelihoodMessage} # up receive message buffer (multiple children, multiple messages)
  downRx::Union{Nothing, LikelihoodMessage} # down receive message buffer (one parent)
  upTx::Union{Nothing, LikelihoodMessage} # RESERVED up outgoing message buffer (one parent)
  downTx::Union{Nothing, LikelihoodMessage} # RESERVED down outgoing message buffer (multiple children but one message)
end
MessageStore() = MessageStore(Dict{Int, LikelihoodMessage}(), nothing, nothing, nothing)

## ^^^ DELETE



"""
$(TYPEDEF)

Data structure for each clique in the Bayes (Junction) tree.
"""
mutable struct BayesTreeNodeData
  frontalIDs::Vector{Symbol}
  separatorIDs::Vector{Symbol}
  inmsgIDs::Vector{Symbol} # Int
  potIDs::Vector{Symbol} # Int # this is likely redundant TODO -- remove
  potentials::Vector{Symbol}
  partialpotential::Vector{Bool}

  dwnPotentials::Vector{Symbol}
  dwnPartialPotential::Vector{Bool}

  cliqAssocMat::Array{Bool,2}
  cliqMsgMat::Array{Bool,2}
  directvarIDs::Vector{Symbol}
  directFrtlMsgIDs::Vector{Symbol}
  msgskipIDs::Vector{Symbol}
  itervarIDs::Vector{Symbol}
  directPriorMsgIDs::Vector{Symbol}
  debug
  debugDwn

  allmarginalized::Bool
  initialized::Symbol
  upsolved::Bool
  downsolved::Bool
  isCliqReused::Bool             # holdover

  # keep the Condition and Channel{Int}'s for now
  solveCondition::Condition

  # FIXME consolidate Dict with LikelihoodMessage, (pull model #674)
  solvableDims::Channel{Dict{Symbol, Float64}}

  # Consolidation for #459 complete!
  upMsgChannel::Channel{LikelihoodMessage}
  dwnMsgChannel::Channel{LikelihoodMessage}

  ## DF THIS MUST BE DELETED AND ONLY USE BTND.up[/dwn]MsgChannel 
  # JT Local messages saved for cache and debuging 
  messages::MessageStore
end

function BayesTreeNodeData(;frontalIDs=Symbol[],
                            separatorIDs=Symbol[],
                            inmsgIDs=Symbol[],
                            potIDs=Symbol[],
                            potentials=Symbol[],
                            partialpotential=Bool[],
                            dwnPotentials=Symbol[],
                            dwnPartialPotential=Bool[],
                            cliqAssocMat=Array{Bool}(undef, 0,0),
                            cliqMsgMat=Array{Bool}(undef, 0,0),
                            directvarIDs=Int[],
                            directFrtlMsgIDs=Int[],
                            msgskipIDs=Int[],
                            itervarIDs=Int[],
                            directPriorMsgIDs=Int[],
                            debug=nothing,
                            debugDwn=nothing,
                            allmarginalized=false,
                            initialized=:null,
                            upsolved=false,
                            downsolved=false,
                            isCliqReused=false,
                            solveCondition=Condition(),
                            solvableDims=Channel{Dict{Symbol,Float64}}(1),
                            upMsgChannel=Channel{LikelihoodMessage}(1),
                            dwnMsgChannel=Channel{LikelihoodMessage}(1),
                            messages = MessageStore()
                          )
  btnd = BayesTreeNodeData(frontalIDs,
                        separatorIDs,
                        inmsgIDs,
                        potIDs,
                        potentials,
                        partialpotential,
                        dwnPotentials,
                        dwnPartialPotential,
                        cliqAssocMat,
                        cliqMsgMat,
                        directvarIDs,
                        directFrtlMsgIDs,
                        msgskipIDs,
                        itervarIDs,
                        directPriorMsgIDs,
                        debug,
                        debugDwn,
                        allmarginalized,
                        initialized,
                        upsolved,
                        downsolved,
                        isCliqReused,
                        solveCondition,
                        solvableDims,
                        upMsgChannel,
                        dwnMsgChannel,
                        messages  )
  #
  put!(btnd.upMsgChannel, LikelihoodMessage())
  return btnd
end
#


function compare( c1::BayesTreeNodeData,
                  c2::BayesTreeNodeData;
                  skip::Vector{Symbol}=[] )
  #
  TP = true

  TP = TP && c1.frontalIDs == c2.frontalIDs
  TP = TP && c1.separatorIDs == c2.separatorIDs
  TP = TP && c1.inmsgIDs == c2.inmsgIDs
  TP = TP && c1.potIDs == c2.potIDs
  TP = TP && c1.potentials == c2.potentials
  TP = TP && c1.partialpotential == c2.partialpotential
  TP = TP && c1.dwnPotentials == c2.dwnPotentials
  TP = TP && c1.dwnPartialPotential == c2.dwnPartialPotential
  TP = TP && c1.cliqAssocMat == c2.cliqAssocMat
  TP = TP && c1.cliqMsgMat == c2.cliqMsgMat
  TP = TP && c1.directvarIDs == c2.directvarIDs
  TP = TP && c1.directFrtlMsgIDs == c2.directFrtlMsgIDs
  TP = TP && c1.msgskipIDs == c2.msgskipIDs
  TP = TP && c1.itervarIDs == c2.itervarIDs
  TP = TP && c1.directPriorMsgIDs == c2.directPriorMsgIDs
  TP = TP && c1.debug == c2.debug
  TP = TP && c1.debugDwn == c2.debugDwn
  TP = TP && c1.allmarginalized == c2.allmarginalized
  TP = TP && c1.initialized == c2.initialized
  TP = TP && c1.upsolved == c2.upsolved
  TP = TP && c1.downsolved == c2.downsolved
  TP = TP && c1.isCliqReused == c2.isCliqReused
  TP = TP && getMsgUpThis(c1) == getMsgUpThis(c2)
  # TP = TP && c1.solveCondition == c2.solveCondition
  TP = TP && c1.solvableDims == c2.solvableDims
  TP = TP && getMsgUpChannel(c1) == getMsgUpChannel(c2)
  TP = TP && c1.dwnMsgChannel == c2.dwnMsgChannel

  return TP
end


## Packed types for serialization


mutable struct PackedBayesTreeNodeData
  frontalIDs::Vector{Symbol}
  separatorIDs::Vector{Symbol}
  inmsgIDs::Vector{Symbol} # Int
  potIDs::Vector{Symbol} # Int # this is likely redundant TODO -- remove
  potentials::Vector{Symbol}
  partialpotential::Vector{Bool}
  dwnPotentials::Vector{Symbol}
  dwnPartialPotential::Vector{Bool}
  cliqAssocMat::Array{Bool,2}
  cliqMsgMat::Array{Bool,2}
  directvarIDs::Vector{Symbol} # Int
  directFrtlMsgIDs::Vector{Symbol} # Int
  msgskipIDs::Vector{Symbol} # Int
  itervarIDs::Vector{Symbol} # Int
  directPriorMsgIDs::Vector{Symbol} # Int
end



function convert(::Type{PackedBayesTreeNodeData}, btnd::BayesTreeNodeData)
  return PackedBayesTreeNodeData(
    btnd.frontalIDs,
    btnd.separatorIDs,
    btnd.inmsgIDs,
    btnd.potIDs,
    btnd.potentials,
    btnd.partialpotential,
    btnd.dwnPotentials,
    btnd.dwnPartialPotential,
    btnd.cliqAssocMat,
    btnd.cliqMsgMat,
    btnd.directvarIDs,
    btnd.directFrtlMsgIDs,
    btnd.msgskipIDs,
    btnd.itervarIDs,
    btnd.directPriorMsgIDs  )
end


function convert(::Type{BayesTreeNodeData}, pbtnd::PackedBayesTreeNodeData)
  btnd = BayesTreeNodeData()
    btnd.frontalIDs = pbtnd.frontalIDs
    btnd.separatorIDs = pbtnd.separatorIDs
    btnd.inmsgIDs = pbtnd.inmsgIDs
    btnd.potIDs = pbtnd.potIDs
    btnd.potentials = pbtnd.potentials
    btnd.partialpotential = pbtnd.partialpotential
    btnd.dwnPotentials = pbtnd.dwnPotentials
    btnd.dwnPartialPotential = pbtnd.dwnPartialPotential
    btnd.cliqAssocMat = pbtnd.cliqAssocMat
    btnd.cliqMsgMat = pbtnd.cliqMsgMat
    btnd.directvarIDs = pbtnd.directvarIDs
    btnd.directFrtlMsgIDs = pbtnd.directFrtlMsgIDs
    btnd.msgskipIDs = pbtnd.msgskipIDs
    btnd.itervarIDs = pbtnd.itervarIDs
    btnd.directPriorMsgIDs = pbtnd.directPriorMsgIDs
  return btnd
end




##==============================================================================
## Must consolideate tree message
##==============================================================================


function messages(btnd::BayesTreeNodeData)
  @warn("btnd.messages will be deprecated")
  btnd.messages
end

messages(clique::TreeClique) = getCliqueData(clique).messages




## ===========================================================================================
## MUST DEPRECATE BELOW
## ===========================================================================================


"""
$(TYPEDEF)
"""
mutable struct FullExploreTreeType{T, T2, T3 <:InMemoryDFGTypes}
  fg::T3
  bt::T2
  cliq::TreeClique
  prnt::T
  sendmsgs::Vector{LikelihoodMessage}
end

const ExploreTreeType{T} = FullExploreTreeType{T, BayesTree}
const ExploreTreeTypeLight{T} = FullExploreTreeType{T, Nothing}


function ExploreTreeType( fgl::G,
                          btl::AbstractBayesTree,
                          vertl::TreeClique,
                          prt::T,
                          msgs::Array{LikelihoodMessage,1} ) where {G <: AbstractDFG, T}
  #
  ExploreTreeType{T}(fgl, btl, vertl, prt, msgs)
end
