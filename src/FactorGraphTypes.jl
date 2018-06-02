import Base: convert
import Base: ==

@compat abstract type InferenceType end
@compat abstract type PackedInferenceType end

@compat abstract type FunctorInferenceType <: Function end

abstract type InferenceVariable end

# been replaced by Functor types, but may be reused for non-numerical cases
@compat abstract type Pairwise <: InferenceType end
@compat abstract type Singleton <: InferenceType end

@compat abstract type FunctorSingleton <: FunctorInferenceType end
# @compat abstract type FunctorPartialSingleton <: FunctorInferenceType end
@compat abstract type FunctorSingletonNH <: FunctorSingleton end

@compat abstract type FunctorPairwise <: FunctorInferenceType end
@compat abstract type FunctorPairwiseMinimize <: FunctorInferenceType end
@compat abstract type FunctorPairwiseNH <: FunctorPairwise end
# @compat abstract type FunctorPairwiseNHMinimize <: FunctorPairwiseMinimize end # TODO

struct ContinuousScalar <: InferenceVariable
  dims::Int
  ContinuousScalar() = new(1)
end
struct ContinuousMultivariate <:InferenceVariable
  dims::Int
  ContinuousMultivariate() = new()
  ContinuousMultivariate(x) = new(x)
end

const FGG = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Array{Graphs.ExVertex,1},Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}}
const FGGdict = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Dict{Int,Graphs.ExVertex},Dict{Int,Array{Graphs.Edge{Graphs.ExVertex},1}}}

# Condensed representation of KernelDensityEstimate, by saving points and bandwidth
mutable struct EasyMessage
  pts::Array{Float64,2}
  bws::Array{Float64,1}
  EasyMessage() = new()
  EasyMessage(a::Array{Float64,2}, b::Array{Float64,1}) = new(a,b)
  EasyMessage(p::BallTreeDensity) = new(getPoints(p), getBW(p)[:,1])
end


mutable struct FactorGraph
  g::FGGdict
  bn
  IDs::Dict{Symbol,Int}
  fIDs::Dict{Symbol,Int}
  id::Int
  nodeIDs::Array{Int,1} # TODO -- ordering seems improved to use adj permutation -- pending merge JuliaArchive/Graphs.jl/#225
  factorIDs::Array{Int,1}
  bnverts::Dict{Int,Graphs.ExVertex} # TODO -- not sure if this is still used, remove
  bnid::Int # TODO -- not sure if this is still used
  dimID::Int
  cg
  cgIDs::Dict{Int,Int} # cgIDs[exvid] = neoid
  sessionname::String
  robotname::String
  registeredModuleFunctions::VoidUnion{Dict{Symbol, Function}}
  reference::VoidUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}
  stateless::Bool
  FactorGraph() = new()
  FactorGraph(
    x1,
    x2,
    x3,
    x4,
    x5,
    x6,
    x7,
    x8,
    x9,
    x10,
    x11,
    x12,
    x13,
    x14,
    x15,
    x16
   ) = new(
    x1,
    x2,
    x3,
    x4,
    x5,
    x6,
    x7,
    x8,
    x9,
    x10,
    x11,
    x12,
    x13,
    x14,
    x15,
    x16,
    false )
end

"""
    $(SIGNATURES)

Construct an empty FactorGraph object with the minimum amount of information / memory populated.
"""
function emptyFactorGraph(;reference::VoidUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}=nothing)
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
                     Dict{Int,Int}(),
                     "",
                     "",
                     Dict{Symbol, Function}(:IncrementalInference=>IncrementalInference.getSample), # TODO likely to be removed
                     reference  ) #evalPotential
    return fg
end

mutable struct VariableNodeData
  initval::Array{Float64,2}
  initstdev::Array{Float64,2}
  val::Array{Float64,2}
  bw::Array{Float64,2}
  BayesNetOutVertIDs::Array{Int,1}
  dimIDs::Array{Int,1}
  dims::Int
  eliminated::Bool
  BayesNetVertID::Int
  separator::Array{Int,1}
  groundtruth::VoidUnion{ Dict{ Tuple{Symbol, Vector{Float64}} } } # not packed yet
  softtype
  initialized::Bool
  VariableNodeData() = new()
  function VariableNodeData(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
    warn("Deprecated use of VariableNodeData(11 param), use 13 parameters instead")
    new(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11, nothing, true) # TODO ensure this is initialized true is working for most cases
  end
  VariableNodeData(x1::Array{Float64,2},
                   x2::Array{Float64,2},
                   x3::Array{Float64,2},
                   x4::Array{Float64,2},
                   x5::Vector{Int},
                   x6::Vector{Int},
                   x7::Int,
                   x8::Bool,
                   x9::Int,
                   x10::Vector{Int},
                   x11::VoidUnion{ Dict{ Tuple{Symbol, Vector{Float64}} } },
                   x12,
                   x13::Bool) =
    new(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13)
end

mutable struct FactorMetadata
  factoruserdata
  variableuserdata::Union{Vector, Tuple}
  FactorMetadata() = new([], [])
  FactorMetadata(x1, x2::Union{Vector,Tuple}) = new(x1, x2)
end

mutable struct GenericWrapParam{T} <: FunctorInferenceType
  usrfnc!::T
  params::Vector{Array{Float64,2}}
  varidx::Int
  particleidx::Int
  measurement::Tuple #Array{Float64,2}
  samplerfnc::Function # TODO -- remove, since no required. Direct multiple dispatch at solve
  specialzDim::Bool
  partial::Bool
  factormetadata::FactorMetadata
  GenericWrapParam{T}() where {T} = new()
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}) where {T} = new(fnc, t, 1,1, (zeros(0,1),) , +, false, false, FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int) where {T} = new(fnc, t, i, j, (zeros(0,1),) , +, false, false, FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function) where {T} = new(fnc, t, i, j, meas, smpl, false, false, FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool) where {T} = new(fnc, t, i, j, meas, smpl, szd, false, FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool, partial::Bool) where {T} = new(fnc, t, i, j, meas, smpl, szd, partial, FactorMetadata())
end

mutable struct FastRootGenericWrapParam{T} <: Function
  p::Vector{Int}
  perturb::Vector{Float64}
  X::Array{Float64,2}
  Y::Vector{Float64}
  xDim::Int
  zDim::Int
  gwp::GenericWrapParam{T}
  FastRootGenericWrapParam{T}(xArr::Array{Float64,2}, zDim::Int, residfnc::GenericWrapParam{T}) where {T} =
      new(collect(1:size(xArr,1)), zeros(zDim), xArr, zeros(size(xArr,1)), size(xArr,1), zDim, residfnc)
end

mutable struct GenericFunctionNodeData{T, S}
  fncargvID::Array{Int,1}
  eliminated::Bool
  potentialused::Bool
  edgeIDs::Array{Int,1}
  frommodule::S #Union{Symbol, AbstractString}
  fnc::T
  GenericFunctionNodeData{T, S}() where {T, S} = new{T,S}()
  GenericFunctionNodeData{T, S}(x1, x2, x3, x4, x5::S, x6::T) where {T, S} = new{T,S}(x1, x2, x3, x4, x5, x6)
  GenericFunctionNodeData(x1, x2, x3, x4, x5::S, x6::T) where {T, S} = new{T,S}(x1, x2, x3, x4, x5, x6)
end


###

function addGraphsVert!(fgl::FactorGraph,
            exvert::Graphs.ExVertex;
            labels::Vector{<:AbstractString}=String[])
  #
  Graphs.add_vertex!(fgl.g, exvert)
end

function getVertNode(fgl::FactorGraph, id::Int; nt::Symbol=:var, bigData::Bool=false)
  return fgl.g.vertices[id] # check equivalence between fgl.v/f[i] and fgl.g.vertices[i]
  # return nt == :var ? fgl.v[id] : fgl.f[id]
end
function getVertNode(fgl::FactorGraph, lbl::Symbol; nt::Symbol=:var, bigData::Bool=false)
  return getVertNode(fgl, (nt == :var ? fgl.IDs[lbl] : fgl.fIDs[lbl]), nt=nt, bigData=bigData)
end
getVertNode{T <: AbstractString}(fgl::FactorGraph, lbl::T; nt::Symbol=:var, bigData::Bool=false) = getVertNode(fgl, Symbol(lbl), nt=nt, bigData=bigData)



# excessive function, needs refactoring
function updateFullVertData!(fgl::FactorGraph,
    nv::Graphs.ExVertex;
    updateMAPest::Bool=false )
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
function graphsOutNeighbors(fgl::FactorGraph, exVertId::Int; ready::Int=1,backendset::Int=1, needdata::Bool=false)
  graphsOutNeighbors(fgl.g, getVert(fgl,exVertId), ready=ready, backendset=backendset, needdata=needdata)
end

function graphsGetEdge(fgl::FactorGraph, id::Int)
  nothing
end

function graphsDeleteVertex!(fgl::FactorGraph, vert::Graphs.ExVertex)
  warn("graphsDeleteVertex! -- not deleting Graphs.jl vertex id=$(vert.index)")
  nothing
end


#
