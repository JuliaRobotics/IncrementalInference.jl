using MetaGraphs

## Cliques

"""
    $(TYPEDEF)
Structure to store clique data
DEV NOTES: To replace TreeClique completely
    $(FIELDS)
"""
mutable struct TreeClique
  index::Int # see issue #540
  label::Symbol #NOTE this is currently a label such as clique 1, # The drawing label is saved in attributes, JT I'm not sure of the current use
  data::Any#BayesTreeNodeData #FIXME There is circular type usage in TreeClique, BayesTreeNodeData, CliqStateMachineContainer https://github.com/JuliaLang/julia/issues/269
  attributes::Dict{String, Any} #The drawing attributes
  #solveInProgress #on a clique level a "solve in progress" might be very handy
end

TreeClique(i::Int, label::Symbol) = TreeClique(i, label, emptyBTNodeData(), Dict{String,Any}())
TreeClique(i::Int, label::AbstractString) = TreeClique(i, Symbol(label))

Graphs.make_vertex(g::AbstractGraph{TreeClique}, label::AbstractString) = TreeClique(num_vertices(g) + 1, String(label))
Graphs.vertex_index(v::TreeClique) = v.index
Graphs.attributes(v::TreeClique, g::AbstractGraph) = v.attributes

#TODO the label field and label atribute is a bit confusing with accessors.
DFG.getLabel(cliq::TreeClique) = cliq.attributes["label"]
function setLabel!(cliq::TreeClique, lbl::String)
  cliq.attributes["label"] = lbl
  lbl
end

## end Cliques

## Bayes Trees

abstract type AbstractBayesTree end

#TODO - see #540 related to indexing and ids

# BayesTree declarations
const BTGdict = GenericIncidenceList{TreeClique,Edge{TreeClique},Array{TreeClique,1},Array{Array{Edge{TreeClique},1},1}}
"""
$(TYPEDEF)

Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::FactorGraph`.
"""
mutable struct BayesTree <: AbstractBayesTree
  bt::BTGdict
  btid::Int
  cliques::Dict{Int,TreeClique}
  frontals::Dict{Symbol,Int}
  #TEMP JT for evaluation, store message channels associated with edges between nodes Int -> edge id.  TODO rather store in graph
  messages::Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{BeliefMessage},Channel{BeliefMessage}}}}
  variableOrder::Vector{Symbol}
  buildTime::Float64
end

BayesTree() = BayesTree(Graphs.inclist(TreeClique,is_directed=true),
                         0,
                         Dict{Int,TreeClique}(),
                         Dict{AbstractString, Int}(),
                         Dict{Int, NamedTuple{(:upMsg, :downMsg),Tuple{Channel{BeliefMessage},Channel{BeliefMessage}}}}(),
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
Data structure for the Bayes (Junction) tree, which is used for inference and constructed from a given `::FactorGraph`.
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
  msgsUp::Vector{BeliefMessage} #TODO towards consolidated messages
  msgsDown::Vector{BeliefMessage}
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
                                   x11::BTND=emptyBTNodeData(),
                                   x13::SimpleLogger=SimpleLogger(Base.stdout);
                                   x4i::Int = x4.index) where {BTND, G <: AbstractDFG}
  #
  CliqStateMachineContainer{BTND, G, typeof(x2), typeof(x3)}(x1,x2,x3,x4,x4i,x5,x6,x7,x8,x9,x10,x10aa,x10aaa,x10b,x11,x13, BeliefMessage[], BeliefMessage[])
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

  # future might concentrate these four fields down to two
  # these should become specialized BeliefMessage type
  upMsg::TempBeliefMsg # Dict{Symbol, BallTreeDensity}
  dwnMsg::TempBeliefMsg # Dict{Symbol, BallTreeDensity}
  upInitMsgs::Dict{Int, TempBeliefMsg}
  downInitMsg::TempBeliefMsg

  allmarginalized::Bool
  initialized::Symbol
  upsolved::Bool
  downsolved::Bool
  initUpChannel::Channel{Symbol}
  initDownChannel::Channel{Symbol}
  solveCondition::Condition
  lockUpStatus::Channel{Int}
  lockDwnStatus::Channel{Int}
  solvableDims::Channel{Dict{Symbol, Float64}}
  statehistory::Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}
  #  iSAM2 style
  isCliqReused::Bool
  BayesTreeNodeData() = new()
  BayesTreeNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],
                                x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],x[20],
                                x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31],
                                Vector{Tuple{DateTime, Int, Function, CliqStateMachineContainer}}(), false  )
end

# TODO -- this should be a constructor
function emptyBTNodeData()
  BayesTreeNodeData(Symbol[],Symbol[],Symbol[],
                    Symbol[],Symbol[],Bool[], # 6
                    Symbol[],Bool[],
                    Array{Bool}(undef, 0,0),
                    Array{Bool}(undef, 0,0),
                    Int[],Int[],             # 10+2
                    Int[],Int[],Int[],       # 13+2
                    nothing, nothing,        # 15+2
                    Dict{Symbol, BallTreeDensity}(),  # :null => AMP.manikde!(zeros(1,1), [1.0;], (:Euclid,))),
                    Dict{Symbol, BallTreeDensity}(),  # :null => AMP.manikde!(zeros(1,1), [1.0;], (:Euclid,))),
                    Dict{Int, TempBeliefMsg}(),
                    TempBeliefMsg(),         # 19+2
                    false, :null,
                    false, false,            # 23+2
                    Channel{Symbol}(1), Channel{Symbol}(1), Condition(), # 26+2
                    Channel{Int}(1), Channel{Int}(1),
                    Channel{Dict{Symbol,Float64}}(1) )
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
  IDvals::Dict{Symbol, EasyMessage}
  keepupmsgs::TempBeliefMsg # Dict{Symbol, BallTreeDensity} # TODO Why separate upMsgs?
  totalsolve::Bool
  UpReturnBPType() = new()
  UpReturnBPType(x1,x2,x3,x4,x5) = new(x1,x2,x3,x4,x5)
end

"""
$(TYPEDEF)

TODO refactor msgs into only a single variable
"""
mutable struct DownReturnBPType
  dwnMsg::NBPMessage
  dbgDwn::DebugCliqMCMC
  IDvals::Dict{Symbol,EasyMessage} # Int
  keepdwnmsgs::TempBeliefMsg # Dict{Symbol, BallTreeDensity}
end

"""
$(TYPEDEF)
"""
mutable struct FullExploreTreeType{T, T2, T3 <:InMemoryDFGTypes}
  fg::T3
  bt::T2
  cliq::TreeClique
  prnt::T
  sendmsgs::Vector{NBPMessage}
end

const ExploreTreeType{T} = FullExploreTreeType{T, BayesTree}
const ExploreTreeTypeLight{T} = FullExploreTreeType{T, Nothing}


function ExploreTreeType(fgl::G,
                         btl::AbstractBayesTree,
                         vertl::TreeClique,
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
  cliq::TreeClique
  vid::Symbol # Int
  msgs::Array{NBPMessage,1}
  N::Int
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
  btnd = emptyBTNodeData()
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
