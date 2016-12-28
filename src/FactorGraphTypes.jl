import Base.convert
import Base.==

abstract InferenceType
abstract PackedInferenceType

abstract Pairwise <: InferenceType
abstract Singleton <: InferenceType

abstract FunctorInferenceType <: Function
abstract FunctorPairwise <: FunctorInferenceType

typealias FGG Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Array{Graphs.ExVertex,1},Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}}
typealias FGGdict Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Dict{Int,Graphs.ExVertex},Dict{Int,Array{Graphs.Edge{Graphs.ExVertex},1}}}


type FactorGraph
  g::FGGdict
  bn
  IDs::Dict{Symbol,Int}
  fIDs::Dict{Symbol,Int}
  id::Int64
  nodeIDs::Array{Int,1} # TODO -- ordering seems improved to use adj permutation -- pending merge JuliaArchive/Graphs.jl/#225
  factorIDs::Array{Int,1}
  bnverts::Dict{Int,Graphs.ExVertex} # TODO -- not sure if this is still used, remove
  bnid::Int # TODO -- not sure if this is still used
  dimID::Int64
  cg
  cgIDs::Dict{Int64,Int64}
  sessionname::AbstractString
  registeredModuleFunctions::Dict{Symbol, Function}
  FactorGraph() = new()
  FactorGraph(x...) = new(
    x[1],
    x[2],
    x[3],
    x[4],
    x[5],
    x[6],
    x[7],
    x[8],
    x[9],
    x[10],
    x[11],
    x[12],
    x[13],
    x[14] )
    # x[3] )
    # x[4] ) # removed fg.v
end

function emptyFactorGraph()
    fg = FactorGraph(Graphs.incdict(Graphs.ExVertex,is_directed=false),
                     Graphs.incdict(Graphs.ExVertex,is_directed=true),
                    #  Dict{Int,Graphs.ExVertex}(),
                    #  Dict{Int,Graphs.ExVertex}(),
                     Dict{Symbol,Int}(),
                     Dict{Symbol,Int}(),
                     0,
                     [],
                     [],
                     Dict{Int,Graphs.ExVertex}(),
                     0,
                     0,
                     nothing,
                     Dict{Int64,Int64}(),
                     "",
                     Dict{Symbol, Function}(:IncrementalInference=>IncrementalInference.evalPotential) )
    return fg
end

type VariableNodeData
  initval::Array{Float64,2}
  initstdev::Array{Float64,2}
  val::Array{Float64,2}
  bw::Array{Float64,2}
  BayesNetOutVertIDs::Array{Int64,1}
  dimIDs::Array{Int64,1}
  dims::Int64
  eliminated::Bool
  BayesNetVertID::Int64
  separator::Array{Int64,1}
  VariableNodeData() = new()
  VariableNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10])
end

type PackedVariableNodeData
  vecinitval::Array{Float64,1}
  diminitval::Int64
  vecinitstdev::Array{Float64,1}
  diminitdev::Int64
  vecval::Array{Float64,1}
  dimval::Int64
  vecbw::Array{Float64,1}
  dimbw::Int64
  BayesNetOutVertIDs::Array{Int64,1}
  dimIDs::Array{Int64,1}
  dims::Int64
  eliminated::Bool
  BayesNetVertID::Int64
  separator::Array{Int64,1}
  PackedVariableNodeData() = new()
  PackedVariableNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14])
end

type GenericFunctionNodeData{T, S}
  fncargvID::Array{Int64,1}
  eliminated::Bool
  potentialused::Bool
  edgeIDs::Array{Int64,1}
  frommodule::S #Union{Symbol, AbstractString}
  fnc::T
  GenericFunctionNodeData() = new()
  GenericFunctionNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6])
end

typealias FunctionNodeData{T <: Union{InferenceType, FunctorInferenceType}} GenericFunctionNodeData{T, Symbol}
FunctionNodeData() = GenericFunctionNodeData{T, Symbol}()
FunctionNodeData(x...) = GenericFunctionNodeData{T, Symbol}(x[1],x[2],x[3],x[4],x[5],x[6])

typealias PackedFunctionNodeData{T <: PackedInferenceType} GenericFunctionNodeData{T, AbstractString}
PackedFunctionNodeData() = GenericFunctionNodeData{T, AbstractString}()
PackedFunctionNodeData(x...) = GenericFunctionNodeData{T, AbstractString}(x[1],x[2],x[3],x[4],x[5],x[6])


function convert(::Type{PackedVariableNodeData}, d::VariableNodeData)
  return PackedVariableNodeData(d.initval[:],size(d.initval,1),
                              d.initstdev[:],size(d.initstdev,1),
                              d.val[:],size(d.val,1),
                              d.bw[:],size(d.bw,1),
                              d.BayesNetOutVertIDs,
                              d.dimIDs, d.dims, d.eliminated,
                              d.BayesNetVertID, d.separator)
end
function convert(::Type{VariableNodeData}, d::PackedVariableNodeData)
  r1 = d.diminitval
  c1 = floor(Int,length(d.vecinitval)/r1)
  M1 = reshape(d.vecinitval,r1,c1)

  r2 = d.diminitdev
  c2 = floor(Int,length(d.vecinitstdev)/r2)
  M2 = reshape(d.vecinitstdev,r2,c2)

  r3 = d.dimval
  c3 = floor(Int,length(d.vecval)/r3)
  M3 = reshape(d.vecval,r3,c3)

  r4 = d.dimbw
  c4 = floor(Int,length(d.vecbw)/r4)
  M4 = reshape(d.vecbw,r4,c4)

  return VariableNodeData(M1,M2,M3,M4, d.BayesNetOutVertIDs,
    d.dimIDs, d.dims, d.eliminated, d.BayesNetVertID, d.separator)
end
function VNDencoder(P::Type{PackedVariableNodeData}, d::VariableNodeData)
  return convert(P, d) #PackedVariableNodeData
end
function VNDdecoder(T::Type{VariableNodeData}, d::PackedVariableNodeData)
  return convert(T, d) #VariableNodeData
end


function compare(a::VariableNodeData,b::VariableNodeData)
    TP = true
    TP = TP && a.initval == b.initval
    TP = TP && a.initstdev == b.initstdev
    TP = TP && a.val == b.val
    TP = TP && a.bw == b.bw
    TP = TP && a.BayesNetOutVertIDs == b.BayesNetOutVertIDs
    TP = TP && a.dimIDs == b.dimIDs
    TP = TP && a.dims == b.dims
    TP = TP && a.eliminated == b.eliminated
    TP = TP && a.BayesNetVertID == b.BayesNetVertID
    TP = TP && a.separator == b.separator
    return TP
end

function ==(a::VariableNodeData,b::VariableNodeData, nt::Symbol=:var)
  return IncrementalInference.compare(a,b)
end

function addGraphsVert!{T <: AbstractString}(fgl::FactorGraph, exvert::Graphs.ExVertex;
    labels::Vector{T}=String[])
  Graphs.add_vertex!(fgl.g, exvert)
end

function getVertNode(fgl::FactorGraph, id::Int64, nt::Symbol=:var)
  return fgl.g.vertices[id] # check equivalence between fgl.v/f[i] and fgl.g.vertices[i]
  # return nt == :var ? fgl.v[id] : fgl.f[id]
end
function getVertNode(fgl::FactorGraph, lbl::Symbol, nt::Symbol=:var)
  return getVertNode(fgl, (nt == :var ? fgl.IDs[lbl] : fgl.fIDs[lbl]) , nt)
end
getVertNode{T <: AbstractString}(fgl::FactorGraph, lbl::T, nt::Symbol=:var) = getVertNode(fgl, Symbol(lbl), nt=nt)



# excessive function, needs refactoring
function updateFullVertData!(fgl::FactorGraph,
    nv::Graphs.ExVertex; updateMAPest=false)

  # not required, since we using reference -- placeholder function CloudGraphs interface
  # getVertNode(fgl, nv.index).attributes["data"] = nv.attributes["data"]
  nothing
end


function makeAddEdge!(fgl::FactorGraph, v1::Graphs.ExVertex, v2::Graphs.ExVertex; saveedgeID::Bool=true)
  edge = Graphs.make_edge(fgl.g, v1, v2)
  Graphs.add_edge!(fgl.g, edge)
  if saveedgeID push!(v2.attributes["data"].edgeIDs,edge.index) end
  edge
end

function graphsOutNeighbors(fgl::FactorGraph, vert::Graphs.ExVertex; ready::Int=1,backendset::Int=1)
  Graphs.out_neighbors(vert, fgl.g)
end
function graphsOutNeighbors(fgl::FactorGraph, exVertId::Int64; ready::Int=1,backendset::Int=1)
  graphsOutNeighbors(fgl.g, getVert(fgl,exVertId))
end

function graphsGetEdge(fgl::FactorGraph, id::Int64)
  nothing
end

function graphsDeleteVertex!(fgl::FactorGraph, vert::Graphs.ExVertex)
  warn("graphsDeleteVertex! -- not deleting Graphs.jl vertex id=$(vert.index)")
  nothing
end

#
