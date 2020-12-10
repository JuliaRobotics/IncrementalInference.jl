# Clique types

##==============================================================================
## BayesTreeNodeData
##==============================================================================

mutable struct MessageBuffer
  upRx::Dict{Int, LikelihoodMessage} # up receive message buffer (multiple children, multiple messages)
  downRx::Union{Nothing, LikelihoodMessage} # down receive message buffer (one parent)
  upTx::Union{Nothing, LikelihoodMessage} # RESERVED up outgoing message buffer (one parent)
  downTx::Union{Nothing, LikelihoodMessage} # RESERVED down outgoing message buffer (multiple children but one message)
end
MessageBuffer() = MessageBuffer(Dict{Int, LikelihoodMessage}(), nothing, nothing, nothing)


"""
$(TYPEDEF)

Data structure for each clique in the Bayes (Junction) tree.
"""
mutable struct BayesTreeNodeData
  status::CliqStatus
  frontalIDs::Vector{Symbol}
  separatorIDs::Vector{Symbol}
  inmsgIDs::Vector{Symbol} # Int
  potIDs::Vector{Symbol} # Int # this is likely redundant TODO -- remove
  potentials::Vector{Symbol}
  partialpotential::Vector{Bool}

  dwnPotentials::Vector{Symbol}
  dwnPartialPotential::Vector{Bool}

  cliqAssocMat::Array{Bool,2}
  cliqMsgMat::Array{Bool,2}
  directvarIDs::Vector{Symbol}
  directFrtlMsgIDs::Vector{Symbol}
  msgskipIDs::Vector{Symbol}
  itervarIDs::Vector{Symbol}
  directPriorMsgIDs::Vector{Symbol}
  debug
  debugDwn

  allmarginalized::Bool
  initialized::Symbol
  upsolved::Bool
  downsolved::Bool
  isCliqReused::Bool             # holdover

  # JT Local messages saved for cache and debuging 
  messages::MessageBuffer
end



function BayesTreeNodeData(;status::CliqStatus=NULL,
                            frontalIDs=Symbol[],
                            separatorIDs=Symbol[],
                            inmsgIDs=Symbol[],
                            potIDs=Symbol[],
                            potentials=Symbol[],
                            partialpotential=Bool[],
                            dwnPotentials=Symbol[],
                            dwnPartialPotential=Bool[],
                            cliqAssocMat=Array{Bool}(undef, 0,0),
                            cliqMsgMat=Array{Bool}(undef, 0,0),
                            directvarIDs=Int[],
                            directFrtlMsgIDs=Int[],
                            msgskipIDs=Int[],
                            itervarIDs=Int[],
                            directPriorMsgIDs=Int[],
                            debug=nothing,
                            debugDwn=nothing,
                            allmarginalized=false,
                            initialized=:NULL,
                            upsolved=false,
                            downsolved=false,
                            isCliqReused=false,
                            messages = MessageBuffer()
                          )
  btnd = BayesTreeNodeData(status,
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
                        messages  )
  #
  return btnd
end
#


function compare( c1::BayesTreeNodeData,
                  c2::BayesTreeNodeData;
                  skip::Vector{Symbol}=[] )
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


## Packed types for serialization


mutable struct PackedBayesTreeNodeData
  frontalIDs::Vector{Symbol}
  separatorIDs::Vector{Symbol}
  inmsgIDs::Vector{Symbol} # Int
  potIDs::Vector{Symbol} # Int # this is likely redundant TODO -- remove
  potentials::Vector{Symbol}
  partialpotential::Vector{Bool}
  dwnPotentials::Vector{Symbol}
  dwnPartialPotential::Vector{Bool}
  cliqAssocMat::Array{Bool,2}
  cliqMsgMat::Array{Bool,2}
  directvarIDs::Vector{Symbol} # Int
  directFrtlMsgIDs::Vector{Symbol} # Int
  msgskipIDs::Vector{Symbol} # Int
  itervarIDs::Vector{Symbol} # Int
  directPriorMsgIDs::Vector{Symbol} # Int
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
    btnd.directPriorMsgIDs  )
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
struct CliqueId{T}
  value::T
end

Base.getindex(cId::CliqueId) = cId.value
Base.show(io::IO, ::MIME"text/plain", x::CliqueId) = print(io, x.value)
Base.show(io::IO, x::CliqueId) = print(io, x.value)

"""
    $(TYPEDEF)
Structure to store clique data
DEV NOTES: To replace TreeClique completely
    $(FIELDS)
"""
mutable struct TreeClique
  "Interger id unique within a tree with userId, robotId, sessionId"
  id::CliqueId{Int64} # not to be confused with the underlying index used by LightGraphs.jl, see issue #540
  "Data as `BayesTreeNodeData`"
  data::BayesTreeNodeData 
  "Drawing attributes"
  attributes::Dict{String, Any} 
  #solveInProgress #on a clique level a "solve in progress" might be very handy
end

getId(c::TreeClique) = c.id

function Base.getproperty(x::TreeClique,f::Symbol)
  if f == :index
    Base.depwarn("`TreeCliqe` field `index` is deprecated, use `id`", :getproperty)
    f = :id
  end
  getfield(x,f)
end

function Base.setproperty!(x::TreeClique, f::Symbol, val)
  if f == :index
    Base.depwarn("`TreeCliqe` field `index` is deprecated, use `id`", :setproperty!)
    f = :id
  end
  return setfield!(x, f, convert(fieldtype(typeof(x), f), val))
end

@deprecate TreeClique(i::Int, label::Union{AbstractString, Symbol}) TreeClique(CliqueId(i), BayesTreeNodeData(), Dict{String,Any}())
TreeClique(i::Int) = TreeClique(CliqueId(i), BayesTreeNodeData(), Dict{String,Any}())
TreeClique(id::CliqueId) = TreeClique(id, BayesTreeNodeData(), Dict{String,Any}())


DFG.getLabel(cliq::TreeClique) = cliq.attributes["label"]
function setLabel!(cliq::TreeClique, lbl::String)
  cliq.attributes["label"] = lbl
  lbl
end



## end Cliques
