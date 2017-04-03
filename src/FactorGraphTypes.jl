import Base.convert
import Base.==

abstract InferenceType
abstract PackedInferenceType

# been replaced by Functor types, but may be reused for non-numerical cases
abstract Pairwise <: InferenceType
abstract Singleton <: InferenceType

abstract FunctorInferenceType <: Function

abstract FunctorSingleton <: FunctorInferenceType
# abstract FunctorPartialSingleton <: FunctorInferenceType
abstract FunctorSingletonNH <: FunctorSingleton

abstract FunctorPairwise <: FunctorInferenceType
abstract FunctorPairwiseMinimize <: FunctorInferenceType
abstract FunctorPairwiseNH <: FunctorPairwise


typealias FGG Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Array{Graphs.ExVertex,1},Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}}
typealias FGGdict Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Dict{Int,Graphs.ExVertex},Dict{Int,Array{Graphs.Edge{Graphs.ExVertex},1}}}

typealias VoidUnion{T} Union{Void, T}

# Condensed representation of KernelDensityEstimate, by saving points and bandwidth
type EasyMessage
  pts::Array{Float64,2}
  bws::Array{Float64,1}
end


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
  cgIDs::Dict{Int64,Int64} # cgIDs[exvid] = neoid
  sessionname::AbstractString
  registeredModuleFunctions::VoidUnion{Dict{Symbol, Function}}
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
                     Dict{Symbol, Function}(:IncrementalInference=>IncrementalInference.getSample) ) #evalPotential
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
  groundtruth::VoidUnion{ Dict{ Tuple{Symbol, Vector{Float64}} } }
  VariableNodeData() = new()
  VariableNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11])
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

type GenericWrapParam{T} <: FunctorInferenceType
  usrfnc!::T
  params::Vector{Array{Float64,2}}
  varidx::Int
  particleidx::Int
  measurement::Tuple #Array{Float64,2}
  samplerfnc::Function # TODO -- remove, since no required. Direct multiple dispatch at solve
  specialzDim::Bool
  partial::Bool
  GenericWrapParam() = new()
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}) = new(fnc, t, 1,1, (zeros(0,1),) , +, false, false)
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int) = new(fnc, t, i, j, (zeros(0,1),) , +, false, false)
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function) = new(fnc, t, i, j, meas, smpl, false, false)
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool) = new(fnc, t, i, j, meas, smpl, szd, false)
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool, partial::Bool) = new(fnc, t, i, j, meas, smpl, szd, partial)
end
type FastRootGenericWrapParam{T} <: Function
  p::Vector{Int}
  perturb::Vector{Float64}
  X::Array{Float64,2}
  Y::Vector{Float64}
  xDim::Int
  zDim::Int
  gwp::GenericWrapParam{T}
  FastRootGenericWrapParam{T}(xArr::Array{Float64,2}, zDim::Int, residfnc::T) =
      new(collect(1:size(xArr,1)), zeros(zDim), xArr, zeros(size(xArr,1)), size(xArr,1), zDim, residfnc)
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
    d.dimIDs, d.dims, d.eliminated, d.BayesNetVertID, d.separator,
    nothing)
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


# heavy use of multiple dispatch for converting between packed and original data types during DB usage
function convert{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  return FunctionNodeData{T}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), convert(T, d.fnc))
end
function convert{P <: PackedInferenceType, T <: InferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return PackedFunctionNodeData{P}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc))
end


# Functor version -- TODO, abstraction can be improved here
function convert{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{GenericWrapParam{T}}}, d::PackedFunctionNodeData{P})
  usrfnc = convert(T, d.fnc)
  gwpf = prepgenericwrapper(Graphs.ExVertex[], usrfnc, getSample)
  return FunctionNodeData{GenericWrapParam{T}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), gwpf) #{T}
end
function convert{P <: PackedInferenceType, T <: FunctorInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return PackedFunctionNodeData{P}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(P, d.fnc.usrfnc!))
end

function FNDencode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
end
function FNDdecode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
end

function FNDencode{T <: InferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
end
function FNDdecode{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
end


# Compare FunctionNodeData
function compare{T,S}(a::GenericFunctionNodeData{T,S},b::GenericFunctionNodeData{T,S})
  # TODO -- beef up this comparison to include the gwp
  TP = true
  TP = TP && a.fncargvID == b.fncargvID
  TP = TP && a.eliminated == b.eliminated
  TP = TP && a.potentialused == b.potentialused
  TP = TP && a.edgeIDs == b.edgeIDs
  TP = TP && a.frommodule == b.frommodule
  TP = TP && typeof(a.fnc) == typeof(b.fnc)
  return TP
end


function addGraphsVert!{T <: AbstractString}(fgl::FactorGraph, exvert::Graphs.ExVertex;
    labels::Vector{T}=String[])
  Graphs.add_vertex!(fgl.g, exvert)
end

function getVertNode(fgl::FactorGraph, id::Int64; nt::Symbol=:var, bigData::Bool=false)
  return fgl.g.vertices[id] # check equivalence between fgl.v/f[i] and fgl.g.vertices[i]
  # return nt == :var ? fgl.v[id] : fgl.f[id]
end
function getVertNode(fgl::FactorGraph, lbl::Symbol; nt::Symbol=:var, bigData::Bool=false)
  return getVertNode(fgl, (nt == :var ? fgl.IDs[lbl] : fgl.fIDs[lbl]), nt=nt , bigData=bigData)
end
getVertNode{T <: AbstractString}(fgl::FactorGraph, lbl::T; nt::Symbol=:var, bigData::Bool=false) = getVertNode(fgl, Symbol(lbl), nt=nt, bigData=bigData)



# excessive function, needs refactoring
function updateFullVertData!(fgl::FactorGraph,
    nv::Graphs.ExVertex;
    updateMAPest::Bool=false)
  #

  # not required, since we using reference -- placeholder function CloudGraphs interface
  # getVertNode(fgl, nv.index).attributes["data"] = nv.attributes["data"]
  nothing
end


function makeAddEdge!(fgl::FactorGraph, v1::Graphs.ExVertex, v2::Graphs.ExVertex; saveedgeID::Bool=true)
  edge = Graphs.make_edge(fgl.g, v1, v2)
  Graphs.add_edge!(fgl.g, edge)
  if saveedgeID push!(getData(v2).edgeIDs,edge.index) end #.attributes["data"]
  edge
end

function graphsOutNeighbors(fgl::FactorGraph, vert::Graphs.ExVertex; ready::Int=1,backendset::Int=1, needdata::Bool=false)
  Graphs.out_neighbors(vert, fgl.g)
end
function graphsOutNeighbors(fgl::FactorGraph, exVertId::Int64; ready::Int=1,backendset::Int=1, needdata::Bool=false)
  graphsOutNeighbors(fgl.g, getVert(fgl,exVertId), ready=ready, backendset=backendset, needdata=needdata)
end

function graphsGetEdge(fgl::FactorGraph, id::Int64)
  nothing
end

function graphsDeleteVertex!(fgl::FactorGraph, vert::Graphs.ExVertex)
  warn("graphsDeleteVertex! -- not deleting Graphs.jl vertex id=$(vert.index)")
  nothing
end


#
