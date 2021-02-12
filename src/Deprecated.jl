
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

# TODO Replace with _CalcFactorParametric API
function (s::MsgPrior{<:ParametricTypes})(X1::AbstractVector{T};
                        userdata::Union{Nothing,FactorMetadata}=nothing) where T<:Real

  if isa(s.Z, Normal)
    meas = s.Z.μ
    σ = s.Z.σ
    #TODO confirm signs
    res = meas - X1[1]
    return (res./σ) .^ 2

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- X1
    return res' * iΣ * res

  else
    #this should not happen
    @error("$s not suported, please use non-parametric")
  end                    #
end

# maybe replace X with a type.
function (s::Prior{<:ParametricTypes})(X1::AbstractVector{T};
                    userdata::Union{Nothing,FactorMetadata}=nothing) where T <: Real

  if isa(s.Z, Normal)
    meas = s.Z.μ
    σ = s.Z.σ
    #TODO confirm signs
    res = meas - X1[1]
    return (res./σ) .^ 2

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- X1
    return res' * iΣ * res # + 2*log(1/(  sqrt(det(Σ)*(2pi)^k) )) ## cancel ×1/2 in calling function ## k = dim(μ)
  else
    #this should not happen
    @error("$s not suported, please use non-parametric")
  end
end


# parametric specific functor
# TODO Replace with _CalcFactorParametric API
function (s::LinearRelative{N,<:ParametricTypes})(
                                X1::AbstractArray{<:Real},
                                X2::AbstractArray{<:Real};
                                userdata::Union{Nothing,FactorMetadata}=nothing ) where N
  #
  # can I change userdata to a keyword arg, DF, No will be resolved with consolidation, #467
  # FIXME, replace if with dispatch
  if isa(s.Z, Normal)
    meas = mean(s.Z)
    σ = std(s.Z)
    # res = similar(X2)
    res = meas - (X2[1] - X1[1])
    return (res/σ) .^ 2

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- (X2 .- X1)
    return res' * iΣ * res

  else
    #this should not happen
    @error("$s not supported, please use non-parametric")
  end
end


# parametric specific functor
# TODO Replace with _CalcFactorParametric API
function (s::EuclidDistance{<:ParametricTypes})(X1::AbstractArray{<:Real},
                                                X2::AbstractArray{<:Real};
                                                userdata::Union{Nothing,FactorMetadata}=nothing )
  #
  # can I change userdata to a keyword arg, DF, No will be resolved with consolidation
  #
  # FIXME, replace if with dispatch
  if isa(s.Z, Normal)
    meas = mean(s.Z)
    σ = std(s.Z)
    # res = similar(X2)
    res = meas - norm(X2 - X1)
    res *= res
    return res/(σ^2)

  elseif isa(s.Z, MvNormal)
    meas = mean(s.Z)
    iΣ = invcov(s.Z)
    #TODO confirm math : Σ^(1/2)*X
    res = meas .- (X2 .- X1)
    return res' * iΣ * res

  else
    #this should not happen
    @error("$s not suported, please use non-parametric")
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

"""
$(TYPEDEF)

TODO TO BE DEPRECATED
"""
mutable struct DebugCliqMCMC
  mcmc::Union{Nothing, Array{CliqGibbsMC,1}}
  outmsg::LikelihoodMessage
  outmsglbls::Dict{Symbol, Symbol} # Int
  priorprods::Vector{CliqGibbsMC}
  DebugCliqMCMC() = new()
  DebugCliqMCMC(a,b,c,d) = new(a,b,c,d)
end




##==============================================================================
## Deprecate code below before v0.22
##==============================================================================



"""
(cf::CalcFactor)( res::AbstractVector{<:Real}, meas..., params...)

Default fallback for the standard factor calculation interface, as in `cf.factor(residual, noise_process, parameters)`,
where factors are either library standard or user out-of-library factor definitions.  See documentation for
more details and tutorials on using your own factors (designed to be as easy as possible).

Notes
- These residual calculations use used to find non-Gaussian / multimodal (incl. PPE) and conventional Gaussian estimates. 
- `cf.legacyMeas == (measparams[1:cf._measCount]...,)`

Example
```julia
# TBD
```
"""
function (cf::CalcFactor)(measparams... ) #where {T<:FunctorInferenceType,M,P<:Tuple,X<:AbstractVector} {T,M,P,X}
  #
  # NOTE this is a legacy interface
  res = zeros(size(cf._legacyMeas,1))
  cf.factor(res, cf.metadata, cf._sampleIdx, cf._legacyMeas, cf._legacyParams...)
end



@deprecate prodmultipleonefullpartials(w...;kw...) prodmultiplefullpartials(w...;kw...)

# """
#     $(SIGNATURES)

# Multiply a single full and several partial dimension constraints.

# DevNotes
# - FIXME consolidate partial and full product AMP API, relates to #1010
# - TODO -- reuse memory rather than rand here
# """
# function prodmultipleonefullpartials( dens::Vector{BallTreeDensity},
#                                       partials::Dict{Int, Vector{BallTreeDensity}},
#                                       Ndims::Int,
#                                       N::Int,
#                                       manis::Tuple  )
#   #

#   # TODO -- should this be [1.0] or ones(Ndims)
#   denspts = getPoints(dens[1])
#   pGM = deepcopy(denspts)

#   for (dimnum,pp) in partials
#     push!(pp, AMP.manikde!(pGM[dimnum:dimnum,:], (manis[dimnum],) ))
#   end
  
#   # do each partial dimension individually
#   for (dimnum,pp) in partials
#     pGM[dimnum,:] = AMP.manifoldProduct(pp, (manis[dimnum],), Niter=1) |> getPoints
#   end
#   return pGM
# end

"""
    $(SIGNATURES)

Multiply different dimensions from partial constraints individually.

DevNotes
- FIXME Integrate with `manifoldProduct`, see #1010
"""
function productpartials!(pGM::Array{Float64,2},
                          dummy::BallTreeDensity,
                          partials::Dict{Int, Vector{BallTreeDensity}},
                          manis::Tuple  )
  #
  @warn "productpartials! is being deprecated without direct replacement, see future versions of AMP.manifoldProduct instead."
  # do each partial dimension individually
  for (dimnum,pp) in partials
    pGM[dimnum,:] = AMP.manifoldProduct(pp, (manis[dimnum],), Niter=1) |> getPoints
  end
  nothing
end

@deprecate freshSamples!(w...;kw...) sampleFactor!(w...;kw...)
@deprecate freshSamples(w...;kw...) sampleFactor(w...;kw...)

@deprecate numericSolutionCCW!(w...;kw...) _solveCCWNumeric!(w...;kw...)

@deprecate solveFactorMeasurements( dfg::AbstractDFG,fctsym::Symbol,solveKey::Symbol=:default;retries::Int=3 ) approxDeconv(dfg,fctsym,solveKey,retries=retries)


# figure out how to deprecate (not critical at the moment)
# used in RoMEPlotting 
# NOTE: TempUpMsgPlotting will be removed.
# Replaced with: UpMsgPlotting = @NamedTuple{cliqId::CliqueId{Int}, depth::Int, belief::TreeBelief}
const TempUpMsgPlotting = Dict{Symbol,Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}}

##==============================================================================
## Deprecate code below before v0.21
##==============================================================================

# see DFG #590
@deprecate extractdistribution(x) convert(SamplableBelief, x)


"""
$(TYPEDEF)

MUST BE DEPRECATED IN FAVOR OF EuclidDistance
"""
mutable struct Ranged <: AbstractRelativeRoots
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    Ranged() = new()
    Ranged(x...) = new(x[1], x[2], x[3])
end
"""
$(TYPEDEF)
"""
mutable struct PackedRanged <: PackedInferenceType
    Zij::Array{Float64,1}
    Cov::Array{Float64,1}
    W::Array{Float64,1}
    PackedRanged() = new()
    PackedRanged(x...) = new(x[1], x[2], x[3])
end
function convert(::Type{Ranged}, r::PackedRanged)
  return Ranged(r.Zij, r.Cov, r.W)
end
function convert(::Type{PackedRanged}, r::Ranged)
  return PackedRanged(r.Zij, r.Cov, r.W)
end
function (ra::Ranged)(res::AbstractVector{<:Real},
    userdata::FactorMetadata,
    idx::Int,
    meas::Tuple{<:AbstractArray{<:Real,2}},
    p1::AbstractArray{<:Real},
    l1::AbstractArray{<:Real})

  res[1] = meas[1][1,idx] - abs(l1[1,idx] - p1[1,idx])
  nothing
end
function getSample(ra::Ranged, N::Int=1)
  @warn("`::Ranged` is being deprecated in favor of `::EuclidDistance`")
  ret = zeros(1,N)
  for i in 1:N
    ret[1,i] = ra.Cov[1]*randn()+ra.Zij[1]
  end
  # rand(Distributions.Normal(odo.Zij[1],odo.Cov[1]), N)'
  return (ret,)
end



function getSample( cfo::CalcFactor, 
                    N::Int=1)
  #
  if !hasfield(typeof(cfo.factor), :specialSampler)
    @warn "`getSample(::MyFactor, ::Int)` API is being deprecated, use `getSample(cf::CalcFactor{<:$(typeof(cfo.factor).name)}, ::Int=1) = cf.factor...` instead.  Similarly for residual factor calculations, see IIF #467."
    getSample(cfo.factor, N)
  else
    @warn "`myfactor.specialSampler` API is being deprecated, use `getSample(cf::CalcFactor{<:$(typeof(cfo.factor).name)}, ::Int=1) = cf.metadata.___` instead.  Similarly for residual factor calculations, see IIF #467."
    cfo.factor.specialSampler(cfo.factor, N, cfo.metadata, cfo.metadata.fullvariables...)
  end
end

function freshSamples(usrfnc::T, 
                      N::Int, 
                      fmd::FactorMetadata, 
                      vnd::Vector=[]) where { T <: FunctorInferenceType }
  #
  @warn("freshSamples(::FunctorInferenceType, ::Int, ::FactorMetadata) has been deprecated, use freshSamples!(::CommonConvWrapper, ::Int) instead.")
  if !hasfield(T, :specialSampler)
    getSample(usrfnc, N)
  else
    usrfnc.specialSampler(usrfnc, N, fmd, vnd...)
  end
end



function calcZDim(usrfnc::T, 
                  Xi::Vector{<:DFGVariable}, 
                  fmd::FactorMetadata=FactorMetadata(Xi, getLabel.(Xi), Vector{Matrix{Float64}}(), :null, nothing) ) where {T <: FunctorInferenceType}
  #
  # # zdim = T != GenericMarginal ? size(getSample(usrfnc, 2)[1],1) : 0
  zdim = if T != GenericMarginal
    error("this version of calcZDim has been made obsolete, use calcZDim(::CalcFactor) instead.")
  #   vnds = Xi # (x->getSolverData(x)).(Xi)
  #   smpls = freshSamples(usrfnc, 2, fmd, vnds)[1]
  #   size(smpls,1)
  else
    0
  end
  return zdim
end

# # TODO make in-place memory version
# function getSample( s::Mixture{N_,F,S,T}, 
#                     N::Int=1,
#                     special...;
#                     kw...  ) where {N_,F<:FunctorInferenceType,S,T}
#   #
#   # TODO consolidate #927, case if mechanics has a special sampler
#   # TODO slight bit of waste in computation, but easiest way to ensure special tricks in s.mechanics::F are included
#   ## example case is old FluxModelsPose2Pose2 requiring velocity
#   smplLambda = hasfield(typeof(s.mechanics), :specialSampler) ? ()->s.specialSampler(s.mechanics, N, special...; kw...)[1] : ()->getSample(s.mechanics, N)[1]
#   smpls = smplLambda()
#     # smpls = Array{Float64,2}(undef,s.dims,N)
#   #out memory should be right size first
#   length(s.labels) != N ? resize!(s, N) : nothing
#   s.labels .= rand(s.diversity, N)
#   for i in 1:N
#     mixComponent = s.components[s.labels[i]]
#     smpls[:,i] = rand(mixComponent,1)
#   end

#   # TODO only does first element of meas::Tuple at this stage, see #1099
#   (smpls, ) # s.labels
# end

# # should not be called in case of Prior
# (s::Mixture)( res::AbstractArray{<:Real},
#               fmd::FactorMetadata,
#               idx::Int,
#               meas::Tuple,
#               X... ) = s.mechanics(res, fmd, idx, meas, X...)
#

# function (s::LinearRelative)( res::AbstractArray{<:Real},
#                               userdata::FactorMetadata,
#                               idx::Int,
#                               meas::Tuple,
#                               X1::AbstractArray{<:Real,2},
#                               X2::AbstractArray{<:Real,2}  )
#   #
#   res[:] = meas[1][:,idx] - (X2[:,idx] - X1[:,idx])
#   nothing
# end



# function (s::EuclidDistance)( res::AbstractArray{<:Real},
#                               fmd::FactorMetadata,
#                               idx::Int,
#                               meas::Tuple,
#                               X1::AbstractArray{<:Real,2},
#                               X2::AbstractArray{<:Real,2}  )
#   #
#   res[1] = meas[1][1,idx] - norm(X2[:,idx] - X1[:,idx])
#   res[1] ^= 2
#   return res[1]
# end

function PriorSphere1(mu::Array{Float64}, cov::Array{Float64,2}, W::Vector{Float64})
  @warn "PriorSphere1(mu,cov,W) is deprecated in favor of PriorSphere1(T(...)) -- use for example PriorSphere1(MvNormal(mu, cov))"
  PriorSphere1(MvNormal(mu[:], cov))
end



@deprecate fetchMsgUpThis(cliq::TreeClique) IIF.getMessageBuffer(cliq).upTx


#
