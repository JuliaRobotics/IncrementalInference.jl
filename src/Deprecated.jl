
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


"""
$(TYPEDEF)

TODO TO BE DEPRECATED
"""
mutable struct PotProd
    Xi::Symbol # Int
    prev::Array{Float64,2}
    product::Array{Float64,2}
    potentials::Array{BallTreeDensity,1}
    potentialfac::Vector{Symbol}
end

"""
$(TYPEDEF)

TODO TO BE DEPRECATED
"""
mutable struct CliqGibbsMC
    prods::Array{PotProd,1}
    lbls::Vector{Symbol}
    CliqGibbsMC() = new()
    CliqGibbsMC(a,b) = new(a,b)
end



##==============================================================================
## Deprecate as part of Manifolds.jl consolidation
##==============================================================================


# Legacy support
getManifolds(::InstanceType{Manifolds.Euclidean{Tuple{N}, ℝ}}) where N = tuple([:Euclid for i in 1:N]...)
getManifolds(::InstanceType{Manifolds.Circle{ℝ}})  = (:Circular,)




##==============================================================================
## Deprecate code below before v0.23
##==============================================================================


# FIXME, much consolidation required here
# Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousEuclid{1}}) = AMP.Euclid
# Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousEuclid{2}}) = AMP.Euclid2
# Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousEuclid{3}}) = AMP.Euclid3
# Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{ContinuousEuclid{4}}) = AMP.Euclid4

# Base.convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{Circular}) = AMP.Circle1

# convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{LinearRelative{1}}) = AMP.Euclid
# convert(::Type{<:ManifoldsBase.Manifold}, ::InstanceType{LinearRelative{2}}) = AMP.Euclid2

export Sphere1

@warn "Deprecating old use of Sphere1, being replaced by Cicular instead"
const Sphere1 = Circular
# @deprecate Sphere1(w...;kw...) Circular(w...;kw...)

@deprecate PriorSphere1(w...;kw...) PriorCircular(w...;kw...)
@deprecate Sphere1Sphere1(w...;kw...) CircularCircular(w...;kw...)
@deprecate PackedPriorSphere1(w...;kw...) PackedPriorCircular(w...;kw...)
@deprecate PackedSphere1Sphere1(w...;kw...) PackedCircularCircular(w...;kw...)


# """
# $(TYPEDEF)

# TODO TO BE DEPRECATED
# """
# mutable struct DebugCliqMCMC
#   mcmc::Union{Nothing, Array{CliqGibbsMC,1}}
#   outmsg::LikelihoodMessage
#   outmsglbls::Dict{Symbol, Symbol} # Int
#   priorprods::Vector{CliqGibbsMC}
#   DebugCliqMCMC() = new()
#   DebugCliqMCMC(a,b,c,d) = new(a,b,c,d)
# end


#
