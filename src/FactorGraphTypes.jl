import Base.convert
import Base.==

abstract Pairwise
abstract Singleton


typealias FGG Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Array{Graphs.ExVertex,1},Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}}

type FactorGraph
  g::FGG
  bn
  v::Dict{Int,Graphs.ExVertex}
  f::Dict{Int,Graphs.ExVertex}
  IDs::Dict{AbstractString,Int}
  fIDs::Dict{AbstractString,Int}
  id::Int64
  nodeIDs::Array{Int,1} # TODO -- ordering seems brittle
  factorIDs::Array{Int,1}
  bnverts::Dict{Int,Graphs.ExVertex} # TODO -- not sure if this is still used
  bnid::Int # TODO -- not sure if this is still used
  dimID::Int64
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
    x[12] )
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
  separator::Array{Int64,1} # Will bring to hard type soon, don't worry about this one just yet
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

function ==(a::VariableNodeData,b::VariableNodeData)
  return compare(a,b)
end

function getVarNode(fgl::FactorGraph, id::Int64)
  return fgl.v[id]
end
function getVarNode(fgl::FactorGraph, lbl::AbstractString)
  return getVarNode(fgl, fgl.IDs[lbl])
end

function addNewVertInGraph!(fgl::FactorGraph, vert::Graphs.ExVertex, id::Int64, lbl::AbstractString)
  vert.attributes = Graphs.AttributeDict() #fg.v[fg.id]
  vert.attributes["label"] = lbl #fg.v[fg.id]
  fgl.v[id] = vert
  fgl.IDs[lbl] = id # fg.id
  nothing
end

#
