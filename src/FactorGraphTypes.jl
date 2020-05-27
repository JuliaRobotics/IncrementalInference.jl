import Base: convert
import Base: ==



# TODO been replaced by Functor types, but may be reused for non-numerical cases
abstract type Pairwise <: InferenceType end
abstract type Singleton <: InferenceType end

# TODO deprecate with standard null hypothesis only
abstract type FunctorSingletonNH <: FunctorSingleton end
abstract type FunctorPairwiseNH <: FunctorPairwise end
# abstract type FunctorPairwiseNHMinimize <: FunctorPairwiseMinimize end # TODO


const FGG = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Array{Graphs.ExVertex,1},Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}}
const FGGdict = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Dict{Int,Graphs.ExVertex},Dict{Int,Array{Graphs.Edge{Graphs.ExVertex},1}}}

const BeliefArray{T} = Union{Array{T,2}, Adjoint{T, Array{T,2}} }

"""
$(TYPEDEF)

Solver parameters for the DistributedFactoGraph.

Dev Notes
- TODO remove NothingUnion
"""
mutable struct SolverParams <: DFG.AbstractParams
  dimID::Int
  # TODO remove NothingUnion
  registeredModuleFunctions::NothingUnion{Dict{Symbol, Function}} # remove from
  reference::NothingUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}
  stateless::Bool
  qfl::Int # Quasi fixed length
  isfixedlag::Bool # true when adhering to qfl window size for solves
  limitfixeddown::Bool # if true, then fixed lag will also not update marginalized during down (default false)
  # new functions
  incremental::Bool
  upsolve::Bool
  downsolve::Bool
  drawtree::Bool
  showtree::Bool
  drawtreerate::Float64
  dbg::Bool
  async::Bool
  limititers::Int
  N::Int
  multiproc::Bool
  logpath::String
  graphinit::Bool
  treeinit::Bool # still experimental with known errors
  algorithms::Vector{Symbol} # list of algorithms to run [:default] is mmisam
  spreadNH::Float64 # experimental, entropy spread adjustment used for both null hypo cases.
  devParams::Dict{Symbol,String}
  SolverParams(;dimID::Int=0,
                registeredModuleFunctions=nothing,
                reference=nothing,
                stateless::Bool=false,
                qfl::Int=99999999999,
                isfixedlag::Bool=false,
                limitfixeddown::Bool=false,
                incremental::Bool=true,
                upsolve::Bool=true,
                downsolve::Bool=true,
                drawtree::Bool=false,
                showtree::Bool=false,
                drawtreerate::Float64=0.5,
                dbg::Bool=false,
                async::Bool=false,
                limititers::Int=500,
                N::Int=100,
                multiproc::Bool=1 < nprocs(),
                logpath::String="/tmp/caesar/$(now())",
                graphinit::Bool=true,
                treeinit::Bool=false,
                algorithms::Vector{Symbol}=[:default],
                spreadNH::Float64=3.0,
                devParams::Dict{Symbol,String}=Dict{Symbol,String}()
              ) = new(dimID,
                      registeredModuleFunctions,
                      reference,
                      stateless,
                      qfl,
                      isfixedlag,
                      limitfixeddown,
                      incremental,
                      upsolve,
                      downsolve,
                      drawtree,
                      showtree,
                      drawtreerate,
                      dbg,
                      async,
                      limititers,
                      N,
                      multiproc,
                      logpath,
                      graphinit,
                      treeinit,
                      algorithms,
                      spreadNH,
                      devParams )
  #
end


"""
    $SIGNATURES

Initialize an empty in-memory DistributedFactorGraph `::DistributedFactorGraph` object.
"""
function initfg(dfg::T=InMemDFGType(solverParams=SolverParams());
                                    sessionname="NA",
                                    robotname="",
                                    username="",
                                    cloudgraph=nothing)::T where T <: AbstractDFG
  #
  return dfg
end


#init an empty fg with a provided type and SolverParams
function initfg(::Type{T}; solverParams=SolverParams(),
                           sessionname="NA",
                           robotname="",
                           username="",
                           cloudgraph=nothing)::AbstractDFG where T <: AbstractDFG
  return T(solverParams=solverParams)
end

function initfg(::Type{T}, solverParams::SolverParams;
                           sessionname="NA",
                           robotname="",
                           username="",
                           cloudgraph=nothing)::AbstractDFG where T <: AbstractDFG
  return T{SolverParams}(solverParams=solverParams)
end

"""
$(TYPEDEF)

TODO remove Union types -- issue #383
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

"""
$(TYPEDEF)
"""
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
  # the actual particle being solved at this moment
  particleidx::Int
  # additional data passed to user function -- optionally used by user function
  factormetadata::FactorMetadata
  # subsection indices to select which params should be used for this hypothesis evaluation
  activehypo::Union{UnitRange{Int},Vector{Int}}
  # a permutation vector for low-dimension solves (FunctorPairwise only)
  p::Vector{Int}
  # slight numerical perturbation for degenerate solver cases such as division by zero
  perturb::Vector{Float64}
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
mutable struct CommonConvWrapper{T} <: FactorOperationalMemory where {T<:FunctorInferenceType}
  ### Values consistent across all threads during approx convolution
  usrfnc!::T # user factor / function
  # general setup
  xDim::Int
  zDim::Int
  # special case settings
  specialzDim::Bool # is there a special zDim requirement -- defined by user
  partial::Bool # is this a partial constraint -- defined by user
  # multi hypothesis settings
  hypotheses::Union{Nothing, Distributions.Categorical} # categorical to select which hypothesis is being considered during convolution operation
  certainhypo::Union{Nothing, Vector{Int}}
  # values specific to one complete convolution operation
  params::Vector{Array{Float64,2}} # parameters passed to each hypothesis evaluation event on user function
  varidx::Int # which index is being solved for in params?
  measurement::Tuple # user defined measurement values for each approxConv operation
  threadmodel::Union{Type{SingleThreaded}, Type{MultiThreaded}}
  ### particular convolution computation values per particle idx (varies by thread)
  cpt::Vector{ConvPerThread}

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



#
