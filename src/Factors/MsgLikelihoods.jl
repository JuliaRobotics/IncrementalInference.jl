
"""
$(TYPEDEF)

Message prior on all dimensions of a variable node in the factor graph.

Notes
- Only temporary existance during CSM operations.
"""
struct MsgPrior{T <: SamplableBelief} <: AbstractPrior
  Z::T
  inferdim::Float64
end

# MsgPrior{T}() where {T} = new{T}()
# MsgPrior{T}(z::T, infd::R) where {T <: SamplableBelief, R <: Real} = new{T}(z, infd)
# function MsgPrior(z::T, infd::R) where {T <: SamplableBelief, R <: Real}
#     MsgPrior{T}(z, infd)
# end

getSample(cf::CalcFactor{<:MsgPrior}, N::Int=1) = (reshape(rand(cf.factor.Z,N),:,N), )


# this is a developmental type, will be standardized after conclusion of #1010
const MsgRelativeType = Vector{NamedTuple{(:variables, :likelihood), Tuple{Vector{Symbol},DFG.AbstractRelative}}}

const MsgPriorType = Dict{Symbol, MsgPrior{BallTreeDensity}}


function (s::CalcFactor{<:MsgPrior})( res::AbstractVector{<:Real},
                                            z,
                                            x1  ) # where {M<:FactorMetadata,P<:Tuple,X<:AbstractVector}
  #
  res .= z .- x1
  nothing
end



struct PackedMsgPrior <: PackedInferenceType
  Z::String
  inferdim::Float64
  # PackedMsgPrior() = new()
  # PackedMsgPrior(z::S, infd::R) where {S <: AbstractString, R <: Real} = new(string(z), infd)
end

function convert(::Type{PackedMsgPrior}, d::MsgPrior)
  PackedMsgPrior(convert(PackedSamplableBelief, d.Z), d.inferdim)
end
function convert(::Type{<:MsgPrior}, d::PackedMsgPrior)
  MsgPrior(convert(SamplableBelief, d.Z), d.inferdim)
end





"""
    $TYPEDEF

Internal development types used during consolidation.  Stores relative and prior information making up a joint likelihood 
message passed upward on the Bayes tree.
"""
mutable struct _MsgJointLikelihood
  relatives::IIF.MsgRelativeType
  priors::IIF.MsgPriorType
end

_MsgJointLikelihood(;relatives::MsgRelativeType=MsgRelativeType(),priors::MsgPriorType=MsgPriorType() ) = _MsgJointLikelihood(relatives, priors)


function Base.show(io::IO, x::_MsgJointLikelihood)
  println(io, )
  printstyled(io, " _MsgJointLikelihood:\n", color=:blue)
  print(io, "  .relatives: ")
  for tp in x.relatives
    print(io, tp.variables, "::", typeof(tp.likelihood).name)
    print(io, "; ")
  end
  println(io)
  print(io, "  .priors: ")
  for k in keys(x.priors)
    print(io, k, ", ")
  end
  println(io)
end

Base.show(io::IO, ::MIME"text/plain", x::_MsgJointLikelihood) = show(io, x)


"""
  $(TYPEDEF)
Belief message for message passing on the tree.  This should be considered an incomplete joint probility.

Notes:
- belief -> Dictionary of [`TreeBelief`](@ref)
- variableOrder -> Ordered variable id list of the seperators in cliqueLikelihood
- cliqueLikelihood -> marginal distribution (<: `SamplableBelief`) over clique seperators.
- Older names include: productFactor, Fnew, MsgPrior, LikelihoodMessage

DevNotes:
- Used by both nonparametric and parametric.
- Objective for parametric case: `MvNormal(μ=[:x0;:x2;:l5], Σ=[+ * *; * + *; * * +])`.
- Part of the consolidation effort, see #459.
- Better conditioning for joint structure in the works using deconvolution, see #579, #635.
  - TODO confirm why <: Singleton.

$(TYPEDFIELDS)
"""
mutable struct LikelihoodMessage{T <: MessageType} <: AbstractPrior
  status::CliqStatus
  belief::Dict{Symbol, TreeBelief} # will eventually be deprecated
  variableOrder::Vector{Symbol}
  cliqueLikelihood::Union{Nothing,SamplableBelief}
  msgType::T
  hasPriors::Bool
  # this is different from belief[].inferdim, as the total available infer dims remaining during down msgs -- see #910
  childSolvDims::Dict{Int, Float64} 
  # calc differential factors for joint in the child clique
  jointmsg::_MsgJointLikelihood
  # diffJoints::Vector{NamedTuple{(:variables, :likelihood), Tuple{Vector{Symbol},DFG.AbstractRelative}}}
end


LikelihoodMessage(; status::CliqStatus=NULL,
                    beliefDict::Dict=Dict{Symbol, TreeBelief}(),
                    variableOrder::Vector{Symbol}=Symbol[],
                    cliqueLikelihood=nothing,
                    msgType::T=NonparametricMessage(),
                    hasPriors::Bool=true,
                    childSolvDims::Dict{Int, Float64}=Dict{Int, Float64}(), 
                    jointmsg::_MsgJointLikelihood=_MsgJointLikelihood(),
                  ) where {T <: MessageType} =
        LikelihoodMessage{T}(status, beliefDict, variableOrder, cliqueLikelihood, msgType, hasPriors, childSolvDims, jointmsg)
#

function Base.show(io::IO, ::MIME"text/plain", msg::LikelihoodMessage)
  t = typeof(msg)
  fields = fieldnames(t)
  nf = nfields(msg)
  println(io, "LikelihoodMessage:")
  for f in fields
    printstyled(io, f,": ", color=:blue)
    show(io, getproperty(msg, f))
    println(io)
  end

end

function compare( l1::LikelihoodMessage,
                  l2::LikelihoodMessage;
                  skip::Vector{Symbol}=[] )
  #
  TP = true
  TP = TP && l1.status == l2.status
  TP = TP && l1.variableOrder == l2.variableOrder
  TP = TP && l1.msgType == l2.msgType
  TP = TP && l1.cliqueLikelihood |> typeof == l2.cliqueLikelihood |> typeof
  for (k,v) in l1.belief
    TP = TP && haskey(l2.belief, k)
    TP = TP && compare(v, l2.belief[k])
  end
end

==(l1::LikelihoodMessage,l2::LikelihoodMessage) = compare(l1,l2)



## =========================================================================================
## DEPRECATE BELOW AS ABLE
## =========================================================================================


# figure out how to deprecate (not critical at the moment)
# used in RoMEPlotting
const TempUpMsgPlotting = Dict{Symbol,Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}}




