


##==============================================================================
## Delete at end v0.12.x
##==============================================================================

function getVertKDE(v::DFGVariable)
  return getKDE(v)
end
function getVertKDE(dfg::AbstractDFG, lbl::Symbol)
  v = DFG.getVariable(dfg, lbl)
  return getKDE(v)
end

@deprecate getVertKDE(v::DFGVariable) getKDE(v)
@deprecate getVertKDE(dfg::AbstractDFG, lbl::Symbol) getKDE(dfg, lbl)


function edgelist2edgedict(edgelist::Array{Graphs.Edge{TreeClique},1})
  error("edgelist2edgedict is obsolete, use DFG methods instead.")
  edgedict = Dict{Int,Graphs.Edge{TreeClique}}()
  for edge in edgelist
    edgedict[edge.index] = edge
  end
  return edgedict
end


# TODO: Confirm this is supposed to be a variable?
function setVal!(v::DFGVariable, em::TreeBelief; solveKey::Symbol=:default)
    @warn "setVal! deprecated, use setValKDE! instead"
    setValKDE!(v, em, solveKey=solveKey)
end
function setVal!(v::DFGVariable, p::BallTreeDensity; solveKey::Symbol=:default)
    @warn "setVal! deprecated, use setValKDE! instead"
    setValKDE!(v, p, solveKey=solveKey)
end

@deprecate manualinit!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity)
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) initManual!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) initManual!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) false



# not sure if and where this is still being used
function _evalType(pt::String)::Type
    @error "_evalType has been deprecated, use DFG serialization methods instead."
    try
        getfield(Main, Symbol(pt))
    catch ex
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
    end
end


##==============================================================================
## Delete at end v0.11.x
##==============================================================================

"""
    $SIGNATURES

writeGraphPdf deprecated, use drawGraph instead
"""
function writeGraphPdf(fgl::AbstractDFG;
                       viewerapp::AbstractString="evince",
                       filepath::AbstractString="/tmp/fg.pdf",
                       engine::AbstractString="neato",
                       show::Bool=true )
  #
  @warn "writeGraphPdf deprecated, use drawGraph instead"
  drawGraph(fgl, viewerapp=viewerapp, filepath=filepath, engine=engine, show=show )
end

# excessive function, needs refactoring
# fgl := srcv
function updateFullVertData!(fgl::AbstractDFG,
                             srcv::DFGNode;
                             updatePPE::Bool=false )
  #
  @warn "Deprecated updateFullVertData!, need alternative likely in DFG.mergeGraphVariableData!"

  sym = Symbol(srcv.label)
  isvar = isVariable(fgl, sym)

  dest = isvar ? DFG.getVariable(fgl, sym) : DFG.getFactor(fgl, sym)
  lvd = getSolverData(dest)
  srcvd = getSolverData(srcv)

  if isvar
    if size(lvd.val) == size(srcvd.val)
      lvd.val .= srcvd.val
    else
      lvd.val = srcvd.val
    end
    lvd.bw[:] = srcvd.bw[:]
    lvd.initialized = srcvd.initialized
    lvd.inferdim = srcvd.inferdim
    setSolvedCount!(lvd, getSolvedCount(srcvd))

    if updatePPE
      # set PPE in dest from values in srcv
      # TODO must work for all keys involved
      # dest := srcv
      updatePPE!(fgl, srcv)
      # getVariablePPEs(dest)[:default] = getVariablePPEs(srcv)[:default]
    end
  else
    # assuming nothing to be done
  end

  nothing
end

@deprecate setData!(v::TreeClique, data) setCliqueData!(v,data)


"""
$(TYPEDEF)

NOTE: Deprecated by DistributedFactorGraphs.
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
  FactorGraph(;reference::NothingUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}=nothing, is_directed::Bool=true ) = new(Graphs.incdict(Graphs.ExVertex,is_directed=false),
                      Graphs.incdict(Graphs.ExVertex,is_directed=is_directed),
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


#
