
function _MsgJointLikelihood(;
  relatives::MsgRelativeType = MsgRelativeType(),
  priors::MsgPriorType = MsgPriorType(),
)
  return _MsgJointLikelihood(relatives, priors)
end

function Base.show(io::IO, x::_MsgJointLikelihood)
  println(io)
  printstyled(io, " _MsgJointLikelihood:\n"; color = :blue)
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
  return println(io)
end

Base.show(io::IO, ::MIME"text/plain", x::_MsgJointLikelihood) = show(io, x)

function LikelihoodMessage(;
  sender::NamedTuple{(:id, :step), Tuple{Int, Int}} = (; id = 0, step = 0),
  status::CliqStatus = NULL,
  beliefDict::Dict = Dict{Symbol, TreeBelief}(),
  variableOrder::Vector{Symbol} = Symbol[],
  cliqueLikelihood = nothing,
  msgType::T = NonparametricMessage(),
  hasPriors::Bool = true,
  childSolvDims::Dict{Int, Float64} = Dict{Int, Float64}(),
  jointmsg::_MsgJointLikelihood = _MsgJointLikelihood(),
) where {T <: MessageType}
  return LikelihoodMessage{T}(
    sender,
    status,
    beliefDict,
    variableOrder,
    cliqueLikelihood,
    msgType,
    hasPriors,
    childSolvDims,
    jointmsg,
  )
end
#

function Base.show(io::IO, msg::LikelihoodMessage)
  t = typeof(msg)
  fields = fieldnames(t)
  nf = nfields(msg)
  println(io, "LikelihoodMessage:")
  for f in fields
    printstyled(io, f, ": "; color = :blue)
    show(io, getproperty(msg, f))
    println(io)
  end
end

Base.show(io::IO, ::MIME"text/plain", msg::LikelihoodMessage) = show(io, msg)

function compare(
  l1::LikelihoodMessage,
  l2::LikelihoodMessage;
  skip::Vector{Symbol} = Symbol[],
)
  #
  TP = true
  TP = TP && l1.status == l2.status
  TP = TP && l1.variableOrder == l2.variableOrder
  TP = TP && l1.msgType == l2.msgType
  TP = TP && l1.cliqueLikelihood |> typeof == l2.cliqueLikelihood |> typeof
  for (k, v) in l1.belief
    TP = TP && haskey(l2.belief, k)
    TP = TP && compare(v, l2.belief[k])
  end
  return TP
end

# overload
==(l1::LikelihoodMessage, l2::LikelihoodMessage) = compare(l1, l2)

function BayesTreeNodeData(;
  status::CliqStatus = NULL,
  frontalIDs = Symbol[],
  separatorIDs = Symbol[],
  inmsgIDs = Symbol[],
  potIDs = Symbol[],
  potentials = Symbol[],
  partialpotential = Bool[],
  dwnPotentials = Symbol[],
  dwnPartialPotential = Bool[],
  cliqAssocMat = Array{Bool}(undef, 0, 0),
  cliqMsgMat = Array{Bool}(undef, 0, 0),
  directvarIDs = Int[],
  directFrtlMsgIDs = Int[],
  msgskipIDs = Int[],
  itervarIDs = Int[],
  directPriorMsgIDs = Int[],
  debug = nothing,
  debugDwn = nothing,
  allmarginalized = false,
  initialized = :NULL,
  upsolved = false,
  downsolved = false,
  isCliqReused = false,
  messages = MessageBuffer(),
)
  btnd = BayesTreeNodeData(
    status,
    frontalIDs,
    separatorIDs,
    inmsgIDs,
    potIDs,
    potentials,
    partialpotential,
    dwnPotentials,
    dwnPartialPotential,
    cliqAssocMat,
    cliqMsgMat,
    directvarIDs,
    directFrtlMsgIDs,
    msgskipIDs,
    itervarIDs,
    directPriorMsgIDs,
    debug,
    debugDwn,
    allmarginalized,
    initialized,
    upsolved,
    downsolved,
    isCliqReused,
    messages,
  )
  #
  return btnd
end
#

function compare(c1::BayesTreeNodeData, c2::BayesTreeNodeData; skip::Vector{Symbol} = [])
  #
  TP = true

  TP = TP && c1.frontalIDs == c2.frontalIDs
  TP = TP && c1.separatorIDs == c2.separatorIDs
  TP = TP && c1.inmsgIDs == c2.inmsgIDs
  TP = TP && c1.potIDs == c2.potIDs
  TP = TP && c1.potentials == c2.potentials
  TP = TP && c1.partialpotential == c2.partialpotential
  TP = TP && c1.dwnPotentials == c2.dwnPotentials
  TP = TP && c1.dwnPartialPotential == c2.dwnPartialPotential
  TP = TP && c1.cliqAssocMat == c2.cliqAssocMat
  TP = TP && c1.cliqMsgMat == c2.cliqMsgMat
  TP = TP && c1.directvarIDs == c2.directvarIDs
  TP = TP && c1.directFrtlMsgIDs == c2.directFrtlMsgIDs
  TP = TP && c1.msgskipIDs == c2.msgskipIDs
  TP = TP && c1.itervarIDs == c2.itervarIDs
  TP = TP && c1.directPriorMsgIDs == c2.directPriorMsgIDs
  TP = TP && c1.debug == c2.debug
  TP = TP && c1.debugDwn == c2.debugDwn
  TP = TP && c1.allmarginalized == c2.allmarginalized
  TP = TP && c1.initialized == c2.initialized
  TP = TP && c1.upsolved == c2.upsolved
  TP = TP && c1.downsolved == c2.downsolved
  TP = TP && c1.isCliqReused == c2.isCliqReused

  return TP
end

function convert(::Type{PackedBayesTreeNodeData}, btnd::BayesTreeNodeData)
  return PackedBayesTreeNodeData(
    btnd.frontalIDs,
    btnd.separatorIDs,
    btnd.inmsgIDs,
    btnd.potIDs,
    btnd.potentials,
    btnd.partialpotential,
    btnd.dwnPotentials,
    btnd.dwnPartialPotential,
    btnd.cliqAssocMat,
    btnd.cliqMsgMat,
    btnd.directvarIDs,
    btnd.directFrtlMsgIDs,
    btnd.msgskipIDs,
    btnd.itervarIDs,
    btnd.directPriorMsgIDs,
  )
end

function convert(::Type{BayesTreeNodeData}, pbtnd::PackedBayesTreeNodeData)
  btnd = BayesTreeNodeData()
  btnd.frontalIDs = pbtnd.frontalIDs
  btnd.separatorIDs = pbtnd.separatorIDs
  btnd.inmsgIDs = pbtnd.inmsgIDs
  btnd.potIDs = pbtnd.potIDs
  btnd.potentials = pbtnd.potentials
  btnd.partialpotential = pbtnd.partialpotential
  btnd.dwnPotentials = pbtnd.dwnPotentials
  btnd.dwnPartialPotential = pbtnd.dwnPartialPotential
  btnd.cliqAssocMat = pbtnd.cliqAssocMat
  btnd.cliqMsgMat = pbtnd.cliqMsgMat
  btnd.directvarIDs = pbtnd.directvarIDs
  btnd.directFrtlMsgIDs = pbtnd.directFrtlMsgIDs
  btnd.msgskipIDs = pbtnd.msgskipIDs
  btnd.itervarIDs = pbtnd.itervarIDs
  btnd.directPriorMsgIDs = pbtnd.directPriorMsgIDs
  return btnd
end

##==============================================================================
## Cliques
## TreeClique
##==============================================================================

Base.getindex(cId::CliqueId) = cId.value
Base.show(io::IO, ::MIME"text/plain", x::CliqueId) = print(io, x.value)
Base.show(io::IO, x::CliqueId) = print(io, x.value)

DFG.getId(c::TreeClique) = c.id

TreeClique(i::Int) = TreeClique(CliqueId(i), BayesTreeNodeData(), Dict{String, Any}())
TreeClique(id::CliqueId) = TreeClique(id, BayesTreeNodeData(), Dict{String, Any}())

DFG.getLabel(cliq::TreeClique) = cliq.attributes["label"]
function setLabel!(cliq::TreeClique, lbl::String)
  cliq.attributes["label"] = lbl
  return lbl
end
