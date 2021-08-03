


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




LikelihoodMessage(; sender::NamedTuple{(:id,:step),Tuple{Int,Int}}=(;id=0, step=0),
                    status::CliqStatus=NULL,
                    beliefDict::Dict=Dict{Symbol, TreeBelief}(),
                    variableOrder::Vector{Symbol}=Symbol[],
                    cliqueLikelihood=nothing,
                    msgType::T=NonparametricMessage(),
                    hasPriors::Bool=true,
                    childSolvDims::Dict{Int, Float64}=Dict{Int, Float64}(), 
                    jointmsg::_MsgJointLikelihood=_MsgJointLikelihood(),
                  ) where {T <: MessageType} =
        LikelihoodMessage{T}(sender, status, beliefDict, variableOrder, cliqueLikelihood, msgType, hasPriors, childSolvDims, jointmsg)
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
                  skip::Vector{Symbol}=Symbol[] )
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
  return TP
end

==(l1::LikelihoodMessage,l2::LikelihoodMessage) = compare(l1,l2)


