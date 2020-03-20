
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


##==============================================================================
## Delete in v0.10.x if possible, but definitely by v0.11
##==============================================================================

@deprecate setData!(v::TreeClique, data) setCliqueData!(v,data)


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

@deprecate manualinit!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, vert::DFGVariable, pX::BallTreeDensity)
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) initManual!(dfg::AbstractDFG, sym::Symbol, pX::BallTreeDensity) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) initManual!(dfg::AbstractDFG, sym::Symbol, usefcts::Vector{Symbol}) false
@deprecate manualinit!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) initManual!(dfg::AbstractDFG, sym::Symbol, pts::Array{Float64,2}) false


#
