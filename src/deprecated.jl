# deprecated functions

# must set either dims or initval for proper initialization
# Add node to graph, given graph struct, labal, init values,
# std dev [TODO -- generalize], particle size and ready flag for concurrency
function addNode!(fg::FactorGraph,
      lbl::Symbol,
      initval::Array{Float64}=zeros(1,1),
      stdev::Array{Float64}=ones(1,1); # this is bad and should be removed TODO
      N::Int=100,
      ready::Int=1,
      labels::Vector{<: AbstractString}=String[],
      api::DataLayerAPI=dlapi,
      uid::Int=-1,
      dims::Int=-1  )
  #
  warn("this addNode! is deprecated, please use FactorGraph01.jl:addNode!(fg::FactorGraph, lbl::Symbol, softtype::Type{T}) instead.")
  currid = fg.id+1
  if uid==-1
    fg.id=currid
  else
    currid = uid
  end
  dims = dims != -1 ? dims : size(initval,1)

  lblstr = string(lbl)
  vert = ExVertex(currid,lblstr)
  addNewVarVertInGraph!(fg, vert, currid, lbl, ready, nothing)
  # dlapi.setupvertgraph!(fg, vert, currid, lbl) #fg.v[currid]
  dodims = fg.dimID+1
  # TODO -- vert should not loose information here
  setDefaultNodeData!(vert, initval, stdev, dodims, N, dims) #fg.v[currid]

  vnlbls = deepcopy(labels)
  push!(vnlbls, fg.sessionname)
  # addvert!(fg, vert, api=api)
  api.addvertex!(fg, vert, labels=vnlbls) #fg.g ##vertr =

  fg.dimID+=dims # rows indicate dimensions, move to last dimension
  push!(fg.nodeIDs, currid)

  return vert
end

function getVert(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) where {T <: AbstractString}
  warn("IncrementalInference.getVert{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVert(fgl, Symbol(lbl), api=api)
end

function getVal(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) where {T <: AbstractString}
  warn("IncrementalInference.getVal{T <: AbstractString}(fgl::FactorGraph, lbl::T; api::DataLayerAPI=dlapi) is deprecated, use lbl::Symbol instead")
  getVal(fgl, Symbol(lbl),api=api)
end





function FNDencode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  warn("FNDencode deprecated, use the convert functions through dispatch instead, PackedFunctionNodeData{P=$(P)}.")
  return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
end
function FNDdecode{T <: FunctorInferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  warn("FNDdecode deprecated, use the convert functions through dispatch instead, FunctionNodeData{T=$(T)}.")
  return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
end

function FNDencode{T <: InferenceType, P <: PackedInferenceType}(::Type{PackedFunctionNodeData{P}}, d::FunctionNodeData{T})
  warn("FNDencode deprecated, use the convert functions through dispatch instead, PackedFunctionNodeData{P=$(P)}.")
  return convert(PackedFunctionNodeData{P}, d) #PackedFunctionNodeData{P}
end
function FNDdecode{T <: InferenceType, P <: PackedInferenceType}(::Type{FunctionNodeData{T}}, d::PackedFunctionNodeData{P})
  warn("FNDdecode deprecated, use the convert functions through dispatch instead, FunctionNodeData{T=$(T)}.")
  return convert(FunctionNodeData{T}, d) #FunctionNodeData{T}
end




# this will likely expand with more internal bells and whistles
# to perform in place memory operations for array values in

# see FastRootGenericWrapParam{T}

# TODO -- part of speed and flexibility refactoring exercise
# mutable struct FastGenericRoot{T} <: Function
#   p::Vector{Int}
#   perturb::Vector{Float64}
#   X::Vector{Float64}
#   Y::Vector{Float64}
#   xDim::Int
#   zDim::Int
#   usrfnc::T
#   FastGenericRoot{T}(xDim::Int, zDim::Int, residfnc::T) where {T} =
#       new(collect(1:xDim), zeros(zDim), zeros(xDim), zeros(xDim), xDim, zDim, residfnc)
# end
# function shuffleXAltD!(fr::FastGenericRoot, X::Vector{Float64})
#   copy!(fr.Y, fr.X)
#   for i in 1:fr.zDim
#     fr.Y[fr.p[i]] = X[i]
#   end
#   nothing
# end
# function (fr::FastGenericRoot)( res::Vector{Float64}, x::Vector{Float64} )
#   shuffleXAltD!(fr, x)  #(fr.Y, x, fr.x0, fr.zDim, fr.p)
#   fr.usrfnc( res, fr.Y )
# end
#
# """
#     $(SIGNATURES)
#
# DEPRECATED!
# Solve free variable x by root finding residual function fgr.usrfnc(x, res)
# randomly shuffle x dimensions if underconstrained by measurement z dimensions
# small random perturbation used to prevent trivial solver cases, div by 0 etc.
# result stored in fgr.Y
# """
# function numericRootGenericRandomizedFnc!(
#       fgr::FastGenericRoot{T};
#       perturb::Float64=1e-5,
#       testshuffle::Bool=false ) where {T}
#   #
#   warn("numericRootGenericRandomizedFnc!(fgr::FastGenericRoot{T}...) deprecated, use numericRootGenericRandomizedFnc!(fgr::FastRootGenericWrapRoot{T}...) instead.")
#   fgr.perturb[1:fgr.zDim] = perturb*randn(fgr.zDim)
#   if fgr.zDim < fgr.xDim || testshuffle
#     shuffle!(fgr.p)
#     r = nlsolve(  fgr,
#                 fgr.X[fgr.p[1:fgr.zDim]] + fgr.perturb # this is x0
#              )
#     # copy!(fgr.X, x0) # should need this line?
#     shuffleXAltD!( fgr, r.zero )
#     # copy!(fgr.X, r.zero)
#   else
#     # @show "direct solve"
#     fgr.Y = ( nlsolve(  fgr.usrfnc, fgr.X + fgr.perturb ) ).zero
#     # copy!(fgr.X, y)
#   end
#   nothing
# end
#
#
# function numericRootGenericRandomizedFnc(
#       residFnc!::Function,
#       zDim::Int,
#       xDim::Int,
#       x0::Vector{Float64};
#       perturb::Float64=1e-5,
#       testshuffle::Bool=false   )
#   #
#   # TODO -- this only start of refactoring for inplace, more to come
#   # xDim = length(x0)
#   fgr = FastGenericRoot{typeof(residFnc!)}(xDim, zDim, residFnc!)
#   shuffle!(fgr.p);
#   fgr.perturb[1:fgr.zDim] = perturb*randn(fgr.zDim)
#   copy!(fgr.X, x0)
#
#   numericRootGenericRandomizedFnc!( fgr, perturb=perturb, testshuffle=testshuffle )
#   fgr.Y
# end



# function evalPotentialSpecific(
#       fnc::T,
#       Xi::Vector{Graphs.ExVertex},
#       gwp::GenericWrapParam{T},
#       solvefor::Int;
#       N::Int=100  ) where {T <: FunctorPairwiseMinimize}
#   #
#   # TODO -- could be constructed and maintained at addFactor! time
#   fr, sfidx, maxlen = prepareFastRootGWP(gwp, Xi, solvefor, N)
#   certainidx, allelements, activehypo, mhidx = assembleHypothesesElements!(fr.gwp.hypotheses, maxlen, sfidx, length(Xi))
#
#   # perform the numeric solutions on the indicated elements
#   computeAcrossHypothesis(T, fr, allelements, activehypo, certainidx, sfidx)
#
#   return fr.gwp.params[gwp.varidx]
# end


# function computeAcrossHypothesis(T::Type{<:FunctorPairwiseMinimize}, fr, allelements, activehypo, certainidx, sfidx)
#   count = 0
#   for (mhidx, vars) in activehypo
#     count += 1
#     # if length(allelements[count]) > 0
#     #   fr.gwp.activehypo = vars
#     #   approxConvMinimizeOnElements!(fr, allelements[count])
#     # end
#     if sfidx in certainidx || mhidx in certainidx # certainidx[count] in vars
#       # standard case mhidx, sfidx = $mhidx, $sfidx
#       fr.gwp.activehypo = vars
#       approxConvOnElements!(fr, allelements[count])
#     elseif mhidx == sfidx
#       # multihypo, do conv case, mhidx == sfidx
#       fr.gwp.activehypo = sort(union([sfidx;], certainidx))
#       approxConvOnElements!(fr, allelements[count])
#     elseif mhidx != sfidx
#       # multihypo, take other value case
#       # sfidx=2, mhidx=3:  2 should take a value from 3
#       # sfidx=3, mhidx=2:  3 should take a value from 2
#       fr.gwp.params[sfidx][:,allelements[count]] = fr.gwp.params[mhidx][:,allelements[count]]
#     else
#       error("computeAcrossHypothesis -- not dealing with multi-hypothesis case correctly")
#     end
#   end
#   nothing
# end
