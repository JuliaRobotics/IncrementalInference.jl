


##==============================================================================
## Delete in v0.12
##==============================================================================


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
## Delete in v0.11
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



##==============================================================================
## Delete in v0.10.x if possible, but definitely by v0.11
##==============================================================================

@deprecate setData!(v::TreeClique, data) setCliqueData!(v,data)


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



#
