import Base: convert
import Base: ==

abstract type InferenceType end
abstract type PackedInferenceType end

abstract type FunctorInferenceType <: Function end

abstract type InferenceVariable end
abstract type ConvolutionObject <: Function end

# been replaced by Functor types, but may be reused for non-numerical cases
abstract type Pairwise <: InferenceType end
abstract type Singleton <: InferenceType end

abstract type FunctorSingleton <: FunctorInferenceType end
# abstract type FunctorPartialSingleton <: FunctorInferenceType end
abstract type FunctorSingletonNH <: FunctorSingleton end

abstract type FunctorPairwise <: FunctorInferenceType end
abstract type FunctorPairwiseMinimize <: FunctorInferenceType end
abstract type FunctorPairwiseNH <: FunctorPairwise end
# abstract type FunctorPairwiseNHMinimize <: FunctorPairwiseMinimize end # TODO


const FGG = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Array{Graphs.ExVertex,1},Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}}
const FGGdict = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Dict{Int,Graphs.ExVertex},Dict{Int,Array{Graphs.Edge{Graphs.ExVertex},1}}}

"""
$(TYPEDEF)

Condensed representation of KernelDensityEstimate, by saving points and bandwidth
"""
mutable struct EasyMessage{T <: Tuple}
  pts::Array{Float64,2}
  bws::Array{Float64,1}
  manifolds::T
  EasyMessage{T}() where {T <: Tuple} = new{T}()
  EasyMessage{T}(a::Array{Float64,2}, b::Array{Float64,1}, manis::T) where {T <: Tuple} = new{T}(a,b, manis)
  EasyMessage{T}(p::BallTreeDensity, manis::T) where {T <: Tuple}  = new{T}(getPoints(p), getBW(p)[:,1], manis)
end
EasyMessage(a::Array{Float64,2}, b::Array{Float64,1}, manis::T) where {T <: Tuple} = EasyMessage{T}(a, b, manis)
EasyMessage(p::BallTreeDensity, manis::T) where {T <: Tuple} = EasyMessage{T}(p, manis)


"""
$(TYPEDEF)
"""
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
  username::String
  registeredModuleFunctions::NothingUnion{Dict{Symbol, Function}}
  reference::NothingUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}
  stateless::Bool
  fifo::Vector{Symbol}
  qfl::Int # Quasi fixed length
  isfixedlag::Bool # true when adhering to qfl window size for solves
  FactorGraph(;reference::NothingUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}=nothing ) = new(Graphs.incdict(Graphs.ExVertex,is_directed=false),
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
                      "",
                      Dict{Symbol, Function}(:IncrementalInference=>IncrementalInference.getSample), # TODO likely to be removed
                      reference,
                      false,
                      Symbol[],
                      0,
                      false  )
end

"""
    $SIGNATURES

Initialize an empty `::FactorGraph` object while initializing `sessionname`, `robotname`, and `cloudgraph`.
"""
function initfg(;sessionname="NA",robotname="",username="",cloudgraph=nothing)
  # fgl = RoME.initfg(sessionname=sessionname)
  fgl = FactorGraph()
  fgl.sessionname = sessionname
  fgl.robotname = robotname
  fgl.username = username
  fgl.cg = cloudgraph
  return fgl
end

"""
    $(SIGNATURES)

Construct an empty FactorGraph object with the minimum amount of information / memory populated.

@warn DEPRECATED, use initfg instead.
"""
function emptyFactorGraph(;reference::NothingUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}=nothing)
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
                     "",
                     Dict{Symbol, Function}(:IncrementalInference=>IncrementalInference.getSample), # TODO likely to be removed
                     reference,
                     false,
                     Symbol[],
                     0,
                     false  )
    #
    @warn "emptyFactorGraph(;reference=) is deprecated, use FactorGraph(;reference=) instead."
    return fg
end

"""
$(TYPEDEF)
"""
mutable struct VariableNodeData
  # initval::Array{Float64,2} # TODO deprecate
  # initstdev::Array{Float64,2} # TODO deprecate
  val::Array{Float64,2}
  bw::Array{Float64,2}
  BayesNetOutVertIDs::Array{Int,1}
  dimIDs::Array{Int,1} # Likely deprecate
  dims::Int
  eliminated::Bool
  BayesNetVertID::Int
  separator::Array{Int,1}
  groundtruth::NothingUnion{ Dict{ Tuple{Symbol, Vector{Float64}} } } # not packed yet
  softtype
  initialized::Bool
  ismargin::Bool
  dontmargin::Bool
  VariableNodeData() = new()
  # function VariableNodeData(x1,x2,x3,x4,x5,x6,x7,x8,x9)
  #   @warn "Deprecated use of VariableNodeData(11 param), use 13 parameters instead"
  #   new(x1,x2,x3,x4,x5,x6,x7,x8,x9, nothing, true, false, false) # TODO ensure this is initialized true is working for most cases
  # end
  VariableNodeData(x1::Array{Float64,2},
                   x2::Array{Float64,2},
                   x3::Vector{Int},
                   x4::Vector{Int},
                   x5::Int,
                   x6::Bool,
                   x7::Int,
                   x8::Vector{Int},
                   x9::NothingUnion{ Dict{ Tuple{Symbol, Vector{Float64}} } },
                   x10,
                   x11::Bool,
                   x12::Bool,
                   x13::Bool) =
    new(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13)
end

"""
$(TYPEDEF)
"""
mutable struct FactorMetadata
  factoruserdata
  variableuserdata::Union{Vector, Tuple}
  variablesmalldata::Union{Vector, Tuple}
  solvefor::Union{Symbol, Nothing}
  variablelist::Union{Nothing, Vector{Symbol}}
  dbg::Bool
  FactorMetadata() = new() # [], []
  FactorMetadata(x1, x2::Union{Vector,Tuple},x3) = new(x1, x2, x3, nothing, nothing, false)
  FactorMetadata(x1, x2::Union{Vector,Tuple},x3,x4::Symbol) = new(x1, x2, x3, x4, nothing, false)
  FactorMetadata(x1, x2::Union{Vector,Tuple},x3,x4::Symbol,x5::Vector{Symbol};dbg::Bool=false) = new(x1, x2, x3, x4, x5, dbg)
end

struct SingleThreaded
end
"""
$(TYPEDEF)
"""
struct MultiThreaded
end

"""
$(TYPEDEF)
"""
mutable struct ConvPerThread
  thrid_::Int
  particleidx::Int # the actual particle being solved at this moment
  factormetadata::FactorMetadata # additional data passed to user function -- optionally used by user function
  activehypo::Union{UnitRange{Int},Vector{Int}} # subsection indices to select which params should be used for this hypothesis evaluation
  p::Vector{Int} # a permutation vector for low-dimension solves (FunctorPairwise only)
  perturb::Vector{Float64} # slight numerical perturbation for degenerate solver cases such as division by zero
  X::Array{Float64,2}
  Y::Vector{Float64}
  res::Vector{Float64}
  ConvPerThread() = new()
end

function ConvPerThread(X::Array{Float64,2},
                       zDim::Int;
                       factormetadata::FactorMetadata=FactorMetadata(),
                       particleidx::Int=1,
                       activehypo= 1:length(params),
                       p=collect(1:size(X,1)),
                       perturb=zeros(zDim),
                       Y=zeros(size(X,1)),
                       res=zeros(zDim)  )
  #
  cpt = ConvPerThread()
  cpt.thrid_ = 0
  cpt.X = X
  cpt.factormetadata = factormetadata
  cpt.particleidx = particleidx
  cpt.activehypo = activehypo
  cpt.p = p
  cpt.perturb = perturb
  cpt.Y = Y
  cpt.res = res
  return cpt
end

"""
$(TYPEDEF)
"""
mutable struct CommonConvWrapper{T} <: ConvolutionObject where {T<:FunctorInferenceType}
  ### Values consistent across all threads during approx convolution
  usrfnc!::T # user factor / function
  # general setup
  xDim::Int
  zDim::Int
  # special case settings
  specialzDim::Bool # is there a special zDim requirement -- defined by user
  partial::Bool # is this a partial constraint -- defined by user
  # multi hypothesis settings
  hypotheses::Union{Nothing, Distributions.Categorical} # categorical to select which hypothesis is being considered during convolugtion operation
  certainhypo::Union{Nothing, Vector{Int}}
  # values specific to one complete convolution operation
  params::Vector{Array{Float64,2}} # parameters passed to each hypothesis evaluation event on user function
  varidx::Int # which index is being solved for in params?
  measurement::Tuple # user defined measurement values for each approxConv operation
  threadmodel::Union{Type{SingleThreaded}, Type{MultiThreaded}}

  ### particular convolution computation values per particle idx (varies by thread)
  cpt::Vector{ConvPerThread}
  # varidx::Int # which index is being solved for in params?
  # factormetadata::FactorMetadata # additional data passed to user function -- optionally used by user function
  # activehypo::Union{UnitRange{Int},Vector{Int}} # subsection indices to select which params should be used for this hypothesis evaluation
  # particleidx::Int # the actual particle being solved at this moment
  # p::Vector{Int} # a permutation vector for low-dimension solves (FunctorPairwise only)
  # perturb::Vector{Float64} # slight numerical perturbation for degenerate solver cases such as division by zero
  # X::Array{Float64,2}
  # Y::Vector{Float64}
  # res::Vector{Float64}

  CommonConvWrapper{T}() where {T<:FunctorInferenceType} = new{T}()
end


function CommonConvWrapper(fnc::T,
                           X::Array{Float64,2},
                           zDim::Int,
                           params::Vector{Array{Float64,2}};
                           factormetadata::FactorMetadata=FactorMetadata(),
                           specialzDim::Bool=false,
                           partial::Bool=false,
                           hypotheses=nothing,
                           certainhypo=nothing,
                           activehypo= 1:length(params),
                           varidx::Int=1,
                           measurement::Tuple=(zeros(0,1),),
                           particleidx::Int=1,
                           p=collect(1:size(X,1)),
                           perturb=zeros(zDim),
                           Y=zeros(size(X,1)),
                           xDim=size(X,1),
                           res=zeros(zDim),
                           threadmodel=MultiThreaded  ) where {T<:FunctorInferenceType}
  #
  ccw = CommonConvWrapper{T}()

  ccw.usrfnc! = fnc
  ccw.xDim = xDim
  ccw.zDim = zDim
  ccw.specialzDim = specialzDim
  ccw.partial = partial
  ccw.hypotheses = hypotheses
  ccw.certainhypo=certainhypo
  ccw.params = params
  ccw.varidx = varidx
  ccw.threadmodel = threadmodel
  ccw.measurement = measurement

  # thread specific elements
  ccw.cpt = Vector{ConvPerThread}(undef, Threads.nthreads())
  for i in 1:Threads.nthreads()
    ccw.cpt[i] = ConvPerThread(X, zDim,
                    factormetadata=factormetadata,
                    particleidx=particleidx,
                    activehypo=activehypo,
                    p=p,
                    perturb=perturb,
                    Y=Y,
                    res=res )
  end

  return ccw
end

"""
$(TYPEDEF)
"""
mutable struct GenericFunctionNodeData{T, S}
  fncargvID::Array{Int,1}
  eliminated::Bool
  potentialused::Bool
  edgeIDs::Array{Int,1}
  frommodule::S #Union{Symbol, AbstractString}
  fnc::T
  multihypo::String # likely to moved when GenericWrapParam is refactored
  certainhypo::Vector{Int}
  GenericFunctionNodeData{T, S}() where {T, S} = new{T,S}()
  GenericFunctionNodeData{T, S}(x1, x2, x3, x4, x5::S, x6::T, x7::String="", x8::Vector{Int}=Int[]) where {T, S} = new{T,S}(x1, x2, x3, x4, x5, x6, x7, x8)
  GenericFunctionNodeData(x1, x2, x3, x4, x5::S, x6::T, x7::String="", x8::Vector{Int}=Int[]) where {T, S} = new{T,S}(x1, x2, x3, x4, x5, x6, x7, x8)
  # GenericFunctionNodeData(x1, x2, x3, x4, x5::S, x6::T, x7::String) where {T, S} = new{T,S}(x1, x2, x3, x4, x5, x6, x7)
end


# where {T <: Union{InferenceType, FunctorInferenceType}}
const FunctionNodeData{T} = GenericFunctionNodeData{T, Symbol}
FunctionNodeData(x1, x2, x3, x4, x5::Symbol, x6::T, x7::String="", x8::Vector{Int}=Int[]) where {T <: Union{FunctorInferenceType, ConvolutionObject}}= GenericFunctionNodeData{T, Symbol}(x1, x2, x3, x4, x5, x6, x7, x8)

# where {T <: PackedInferenceType}
const PackedFunctionNodeData{T} = GenericFunctionNodeData{T, <: AbstractString}
PackedFunctionNodeData(x1, x2, x3, x4, x5::S, x6::T, x7::String="", x8::Vector{Int}=Int[]) where {T <: PackedInferenceType, S <: AbstractString} = GenericFunctionNodeData(x1, x2, x3, x4, x5, x6, x7, x8)






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
getVertNode(fgl::FactorGraph, lbl::T; nt::Symbol=:var, bigData::Bool=false) where {T <: AbstractString} = getVertNode(fgl, Symbol(lbl), nt=nt, bigData=bigData)



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
  @warn "graphsDeleteVertex! -- not deleting Graphs.jl vertex id=$(vert.index)"
  nothing
end

function convert(::Type{BallTreeDensity}, p::EasyMessage)
  AMP.manikde!(p.pts, p.bws, p.manifolds)
end

function convert(::Type{EasyMessage}, p::BallTreeDensity, manifolds::T) where {T <: Tuple}
  EasyMessage(getPoints(p), getBW(p)[:,1], manifolds)
end

#
