

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
mutable struct CliqStateMachineContainer{BTND, T <: AbstractDFG, InMemG <: InMemoryDFGTypes, BT <: AbstractBayesTree}
  dfg::T
  cliqSubFg::InMemG
  tree::BT
  cliq::TreeClique
  cliqKey::Int
  parentCliq::Vector{TreeClique}
  childCliqs::Vector{TreeClique}
  forceproceed::Bool # TODO: bad flag that must be removed by refactoring sm
  incremental::Bool
  drawtree::Bool
  dodownsolve::Bool
  delay::Bool
  opts::SolverParams
  refactoring::Dict{Symbol, String}
  oldcliqdata::BTND
  logger::SimpleLogger
  #TODO towards consolidated messages
  # Decision is pull/fetch-model #674 -- i.e. CSM only works inside its own csmc and fetchs messages from neighbors
  # NOTE, DF going for Dict over Vector
  # msgsUp::Dict{Int, LikelihoodMessage} # Vector{LikelihoodMessage}
  # msgsDown::LikelihoodMessage # Vector{LikelihoodMessage}
end

const CSMHistory = Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}

function CliqStateMachineContainer(x1::G,
                                   x2::InMemoryDFGTypes,
                                   x3::AbstractBayesTree,
                                   x4::TreeClique,
                                   x5::Vector{TreeClique},
                                   x6::Vector{TreeClique},
                                   x7::Bool,
                                   x8::Bool,
                                   x9::Bool,
                                   x10::Bool,
                                   x10aa::Bool,
                                   x10aaa::SolverParams,
                                   x10b::Dict{Symbol,String}=Dict{Symbol,String}(),
                                   x11::BTND=BayesTreeNodeData(),
                                   x13::SimpleLogger=SimpleLogger(Base.stdout);
                                   x4i::Int = x4.index) where {BTND, G <: AbstractDFG}
  #
  CliqStateMachineContainer{BTND, G, typeof(x2), typeof(x3)}(x1,x2,x3,x4,x4i,x5,x6,x7,x8,x9,x10,x10aa,x10aaa,x10b,x11,x13 )
              # Dict{Int,LikelihoodMessage}(), LikelihoodMessage())
end



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
  directvarIDs::Vector{Symbol} # Int
  directFrtlMsgIDs::Vector{Symbol} # Int
  msgskipIDs::Vector{Symbol} # Int
  itervarIDs::Vector{Symbol} # Int
  directPriorMsgIDs::Vector{Symbol} # Int
  debug
  debugDwn

  allmarginalized::Bool
  initialized::Symbol
  upsolved::Bool
  downsolved::Bool
  #  iSAM2 style
  isCliqReused::Bool

  # future might concentrate these four fields down to two
  # these should become specialized LikelihoodMessage type
  # TODO, likely to be replaced by Channel counterparts
  upMsg::LikelihoodMessage
  dwnMsg::LikelihoodMessage

  # FIXME Deprecate separate init message locations -- only use up and dwn
  # FIXME ensure these are converted to pull model first #674
  upInitMsgs::Dict{Int, LikelihoodMessage} # FIXME drop dict
  downInitMsg::LikelihoodMessage
  initUpChannel::Channel{LikelihoodMessage}
  initDownChannel::Channel{LikelihoodMessage}

  solveCondition::Condition
  lockUpStatus::Channel{Int}
  lockDwnStatus::Channel{Int}
  # FIXME consolidate Dict with LikelihoodMessage, make pull model first #674
  solvableDims::Channel{Dict{Symbol, Float64}}

  # in and out message channels relating to THIS clique -- only for pull model #674
  upMsgChannel::Channel{LikelihoodMessage}
  dwnMsgChannel::Channel{LikelihoodMessage}
end


function BayesTreeNodeData(;frontalIDs=Symbol[],
                            separatorIDs=Symbol[],
                            inmsgIDs=Symbol[],
                            potIDs=Symbol[],
                            potentials=Symbol[],
                            partialpotential=Bool[],            # 6
                            dwnPotentials=Symbol[],
                            dwnPartialPotential=Bool[],
                            cliqAssocMat=Array{Bool}(undef, 0,0),
                            cliqMsgMat=Array{Bool}(undef, 0,0),
                            directvarIDs=Int[],
                            directFrtlMsgIDs=Int[],             # 10+2
                            msgskipIDs=Int[],
                            itervarIDs=Int[],
                            directPriorMsgIDs=Int[],            # 13+2
                            debug=nothing,
                            debugDwn=nothing,                   # 15+2
                            allmarginalized=false,
                            initialized=:null,
                            upsolved=false,
                            downsolved=false,                   #
                            isCliqReused=false,
                            upMsg=LikelihoodMessage(),
                            dwnMsg=LikelihoodMessage(),
                            upInitMsgs=Dict{Int, LikelihoodMessage}(),
                            downInitMsg=LikelihoodMessage(),         #
                            initUpChannel=Channel{LikelihoodMessage}(1),
                            initDownChannel=Channel{LikelihoodMessage}(1),
                            solveCondition=Condition(),
                            lockUpStatus=Channel{Int}(1),
                            lockDwnStatus=Channel{Int}(1),
                            solvableDims=Channel{Dict{Symbol,Float64}}(1),
                            upMsgChannel=Channel{LikelihoodMessage}(1),
                            dwnMsgChannel=Channel{LikelihoodMessage}(1)
                          )
   BayesTreeNodeData(frontalIDs,
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
                        upMsg,
                        dwnMsg,
                        upInitMsgs,
                        downInitMsg,
                        initUpChannel,
                        initDownChannel,
                        solveCondition,
                        lockUpStatus,
                        lockDwnStatus,
                        solvableDims,
                        upMsgChannel,
                        dwnMsgChannel  )
end
#


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


function ExploreTreeType(fgl::G,
                         btl::AbstractBayesTree,
                         vertl::TreeClique,
                         prt::T,
                         msgs::Array{LikelihoodMessage,1} ) where {G <: AbstractDFG, T}
  #
  ExploreTreeType{T}(fgl, btl, vertl, prt, msgs)
end



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
