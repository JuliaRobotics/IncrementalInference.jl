
##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

export _evalType

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
## TODO deprecated  
##==============================================================================



##==============================================================================
## Deprecate at v0.19
##==============================================================================


@deprecate initManual!(dfg::AbstractDFG, variable::DFGVariable, ptsArr::BallTreeDensity) initManual!(variable, ptsArr)

@deprecate setCliqueDrawColor!(w...;kw...) setCliqueDrawColor!(w...;kw...)

@deprecate evalFactor2(w...;kw...) evalFactor(w...;kw...)


##==============================================================================
## Deprecate at v0.18 
##==============================================================================


# Keep these a bit longer

@deprecate wipeBuildNewTree!(dfg::AbstractDFG; kwargs...) resetBuildTree!(dfg; kwargs...)

@deprecate LinearConditional(N::Int=1) LinearRelative{N}(LinearAlgebra.I)

@deprecate LinearConditional(x...) LinearRelative(x...)

@deprecate PackedLinearConditional(x...) PackedLinearRelative(x...)

@deprecate extractdistribution(x) convert(SamplableBelief, x)


#
