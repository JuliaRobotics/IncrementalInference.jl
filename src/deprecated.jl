# deprecated functions


mutable struct GenericWrapParam{T} <: FunctorInferenceType
  #TODO: <: FunctorIT cannot be right here -- Used in one of the unpacking converters
  usrfnc!::T
  params::Vector{Array{Float64,2}}
  varidx::Int
  particleidx::Int
  measurement::Tuple #Array{Float64,2}
  samplerfnc::Function # TODO -- remove, since no required. Direct multiple dispatch at solve
  specialzDim::Bool
  partial::Bool
  # hypoverts::Vector{Symbol}
  hypotheses::Union{Nothing, Distributions.Categorical}
  activehypo::Union{UnitRange{Int},Vector{Int}}
  factormetadata::FactorMetadata
  GenericWrapParam{T}() where {T} = new()
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}) where {T} = new(fnc, t, 1,1, (zeros(0,1),) , +, false, false, nothing, 1:length(t), FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, varidx::Int, prtcl::Int) where {T} = new(fnc, t, varidx, prtcl, (zeros(0,1),) , +, false, false, nothing, 1:length(t), FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function) where {T} = new(fnc, t, i, j, meas, smpl, false, false, nothing, 1:length(t), FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool) where {T} = new(fnc, t, i, j, meas, smpl, szd, false, nothing, 1:length(t), FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool, partial::Bool) where {T} = new(fnc, t, i, j, meas, smpl, szd, partial, nothing, 1:length(t), FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool, partial::Bool, mhcat::Union{Nothing,Categorical}) where {T} = new(fnc, t, i, j, meas, smpl, szd, partial, mhcat, 1:length(t), FactorMetadata())
  GenericWrapParam{T}(fnc::T, t::Vector{Array{Float64,2}}, i::Int, j::Int, meas::Tuple, smpl::Function, szd::Bool, partial::Bool, mhcat::Tuple) where {T} = new(fnc, t, i, j, meas, smpl, szd, partial, Categorical(mhcat), 1:length(t), FactorMetadata())
end

mutable struct FastRootGenericWrapParam{T} <: Function
  p::Vector{Int}
  perturb::Vector{Float64}
  X::Array{Float64,2}
  Y::Vector{Float64}
  xDim::Int
  zDim::Int
  gwp::GenericWrapParam{T}
  gg::Function
  res::Vector{Float64}
  FastRootGenericWrapParam{T}(xArr::Array{Float64,2}, zDim::Int, residfnc::GenericWrapParam{T}) where {T} =
      new(collect(1:size(xArr,1)), zeros(zDim), xArr, zeros(size(xArr,1)), size(xArr,1), zDim, residfnc, +, zeros(0))
end

function packmultihypo(fnc::GenericWrapParam{T}) where {T<:FunctorInferenceType}
  fnc.hypotheses != nothing ? string(fnc.hypotheses) : ""
end

function prepgenericwrapper(
            Xi::Vector{Graphs.ExVertex},
            usrfnc::UnionAll,
            samplefnc::Function;
            multihypo::Union{Nothing, Distributions.Categorical}=nothing )
      # multiverts::Vector{Symbol}=Symbol[]
  #
  error("prepgenericwrapper -- unknown type usrfnc=$(usrfnc), maybe the wrong usrfnc conversion was dispatched.  Place an error in your unpacking convert function to ensure that IncrementalInference.jl is calling the right unpacking conversion function.")
end

function prepgenericwrapper(
            Xi::Vector{Graphs.ExVertex},
            usrfnc::T,
            samplefnc::Function;
            multihypo::Union{Nothing, Distributions.Categorical}=nothing ) where {T <: FunctorInferenceType}
      # multiverts::Vector{Symbol}=Symbol[]
  #
  warn("GenericWrapParam deprecated, use CommonConvWrapper instead.")
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, 0, 0)
  # test if specific zDim or partial constraint used
  fldnms = fieldnames(usrfnc)
  # sum(fldnms .== :zDim) >= 1
  gwp = GenericWrapParam{T}(
            usrfnc,
            ARR,
            1,
            1,
            (zeros(0,1),),
            samplefnc,
            sum(fldnms .== :zDim) >= 1,
            sum(fldnms .== :partial) >= 1,
            multihypo
        )
    gwp.factormetadata.variableuserdata = []
    gwp.factormetadata.solvefor = :null
    for xi in Xi
      push!(gwp.factormetadata.variableuserdata, getData(xi).softtype)
    end
    return gwp
end

# will be deprecated
function convert(
            ::Type{IncrementalInference.GenericFunctionNodeData{IncrementalInference.GenericWrapParam{F},Symbol}},
            d::IncrementalInference.GenericFunctionNodeData{P,String} ) where {F <: FunctorInferenceType, P <: PackedInferenceType}
  #
  warn("Packing and unpacking of GenericWrapParam{T} will be deprecated, use CommonConvWrapper{T} instead.")
  usrfnc = convert(F, d.fnc)
  # @show d.multihypo
  mhcat = parsemultihypostr(d.multihypo)
  # @show typeof(mhcat)
  gwpf = prepgenericwrapper(Graphs.ExVertex[], usrfnc, getSample, multihypo=mhcat)

  # ccw = prepgenericconvolution(Graphs.ExVertex[], usrfnc, multihypo=mhcat)
  return FunctionNodeData{GenericWrapParam{typeof(usrfnc)}}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), gwpf)
end

function approxConvOnElements!(frl::FastRootGenericWrapParam{T},
                               elements::Union{Vector{Int}, UnitRange{Int}}  ) where {T <: FunctorPairwise}
  #


  for n in elements
    frl.gwp.particleidx = n
    numericRootGenericRandomizedFnc!( frl )
  end
  nothing
end

function approxConvOnElements!(frl::FastRootGenericWrapParam{T},
                               elements::Union{Vector{Int}, UnitRange{Int}}) where {T <: FunctorPairwiseMinimize}

  # TODO should not be claiming new memory every single time....
  # frl.res = zeros(frl.xDim)
  # frl.gg = (x) -> frl.gwp(frl.res, x) # TODO standardize this function into frl


  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  for n in elements
    frl.gwp.particleidx = n
    numericRootGenericRandomizedFnc!( frl )
    # frl.res[1:frl.xDim] = 0.0
    # r = optimize( frl.gg, frl.X[1:frl.xDim, frl.gwp.particleidx] )
    # # TODO -- clearly lots of optmization to be done here
    # frl.Y[1:frl.xDim] = r.minimizer
    # frl.X[1:frl.xDim,frl.gwp.particleidx] = frl.Y
  end
  nothing
end

function prepareFastRootGWP(gwp::GenericWrapParam{T},
                            Xi::Vector{Graphs.ExVertex},
                            solvefor::Int,
                            N::Int  ) where {T <: FunctorInferenceType}
  #


  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
  # should be selecting for the correct multihypothesis mode here with `gwp.params=ARR[??]`
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end
  # Construct complete fr (with fr.gwp) object
  # TODO -- create FastRootGenericWrapParam at addFactor time only?
  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)
  fr.res = zeros(fr.xDim)
  fr.gg = (x) -> fr.gwp(fr.res, x)
  return fr, sfidx, maxlen
end

function computeAcrossHypothesis(frl::FastRootGenericWrapParam{T},
                                 allelements,
                                 activehypo,
                                 certainidx,
                                 sfidx) where {T <:Union{FunctorPairwise, FunctorPairwiseMinimize}}
  #


  count = 0
  for (mhidx, vars) in activehypo
    count += 1
    if sfidx in certainidx || mhidx in certainidx # certainidx[count] in vars
      # standard case mhidx, sfidx = $mhidx, $sfidx
      frl.gwp.activehypo = vars
      approxConvOnElements!(frl, allelements[count])
    elseif mhidx == sfidx
      # multihypo, do conv case, mhidx == sfidx
      frl.gwp.activehypo = sort(union([sfidx;], certainidx))
      approxConvOnElements!(frl, allelements[count])
    elseif mhidx != sfidx
      # multihypo, take other value case
      # sfidx=2, mhidx=3:  2 should take a value from 3
      # sfidx=3, mhidx=2:  3 should take a value from 2
      frl.gwp.params[sfidx][:,allelements[count]] = view(frl.gwp.params[mhidx],:,allelements[count])
      # frl.gwp.params[sfidx][:,allelements[count]] = frl.gwp.params[mhidx][:,allelements[count]]
    else
      error("computeAcrossHypothesis -- not dealing with multi-hypothesis case correctly")
    end
  end
  nothing
end

function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               gwp::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=100,
                               dbg::Bool=false ) where {T <: Union{FunctorPairwise, FunctorPairwiseMinimize}}
  #
  fnc = gwp.usrfnc!

  # Prep computation variables
  fr, sfidx, maxlen = prepareFastRootGWP(gwp, Xi, solvefor, N)
  certainidx, allelements, activehypo, mhidx = assembleHypothesesElements!(fr.gwp.hypotheses, maxlen, sfidx, length(Xi))

  # perform the numeric solutions on the indicated elements
  computeAcrossHypothesis(fr, allelements, activehypo, certainidx, sfidx)

  return fr.gwp.params[gwp.varidx]
end

function assembleNullHypothesis(fr::FastRootGenericWrapParam{T},
                                maxlen::Int,
                                spreadfactor::Float64 ) where {T}
  #


  nhc = rand(fr.gwp.usrfnc!.nullhypothesis, maxlen) - 1
  val = fr.gwp.params[fr.gwp.varidx]
  d = size(val,1)
  var = Base.var(val,2) + 1e-3
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*diagm(var[:]))
  allelements = 1:maxlen
  return allelements, nhc, ENT
end

function computeAcrossNullHypothesis!(frl::FastRootGenericWrapParam{T},
                                      allelements,
                                      nhc,
                                      ENT  ) where {T <: FunctorPairwiseNH}
  #


  # TODO --  Threads.@threads see area4 branch
  for n in allelements
    # frl.gwp(x, res)
    if nhc[n] != 0
      frl.gwp.particleidx = n
      numericRootGenericRandomizedFnc!( frl )
    else
      frl.gwp.params[frl.gwp.varidx][:,n] += rand(ENT)
    end
  end
  nothing
end

function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               generalwrapper::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=100,
                               spreadfactor::Float64=10.0,
                               dbg::Bool=false ) where {T <: FunctorSingletonNH}
  #


  fnc = generalwrapper.usrfnc!

  val = getVal(Xi[1])
  d = size(val,1)
  var = Base.var(val,2) + 1e-3

  # determine amount share of null hypothesis particles
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, N)
  # values of 0 imply null hypothesis
  # generalwrapper.usrfnc!.nullhypothesis::Distributions.Categorical
  nhc = rand(generalwrapper.usrfnc!.nullhypothesis, N) - 1

  # TODO -- not valid for manifold
  ENT = Distributions.MvNormal(zeros(d), spreadfactor*diagm(var[:]))

  for i in 1:N
    if nhc[i] == 0
      generalwrapper.measurement[1][:,i] = val[:,i] + rand(ENT)  # TODO use view and inplace add operation
    end
  end
  # TODO -- returning to memory location inside
  return generalwrapper.measurement[1]
end

function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               generalwrapper::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=0,
                               dbg::Bool=false ) where {T <: FunctorSingleton}
  #


  fnc = generalwrapper.usrfnc!

  nn = (N <= 0 ? size(getVal(Xi[1]),2) : N)
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, nn)
  if !generalwrapper.partial
    return generalwrapper.measurement[1]
  else
    val = deepcopy(getVal(Xi[1]))
    i = 0
    for dimnum in fnc.partial
      i += 1
      val[dimnum,:] = generalwrapper.measurement[1][i,:]
    end
    return val
  end
end

function evalPotentialSpecific(Xi::Vector{Graphs.ExVertex},
                               gwp::GenericWrapParam{T},
                               solvefor::Int;
                               N::Int=100,
                               spreadfactor::Float64=10.0,
                               dbg::Bool=false ) where {T <: FunctorPairwiseNH}
  #

  # TODO -- could be constructed and maintained at addFactor! time
  fr, sfidx, maxlen = prepareFastRootGWP(gwp, Xi, solvefor, N)
  # prepare nullhypothesis
  allelements, nhc, ENT = assembleNullHypothesis(fr, maxlen, spreadfactor)

  # Compute across the true or null hypothesis
  computeAcrossNullHypothesis!(fr, allelements, nhc, ENT )

  return fr.gwp.params[gwp.varidx]
end

# Shuffle incoming X into random permutation in fr.Y
# shuffled fr.Y will be placed back into fr.X[:,fr.gwp.particleidx] upon fr.gwp.usrfnc(x, res)
function shuffleXAltD!(fr::FastRootGenericWrapParam, X::Vector{Float64})

  # populate defaults
  for i in 1:fr.xDim
    fr.Y[i] = fr.X[i, fr.gwp.particleidx]
  end
  # fr.Y[1:fr.xDim] = view(fr.X, 1:fr.xDim, fr.gwp.particleidx)
  # fr.Y[1:fr.xDim] = fr.X[1:fr.xDim,fr.gwp.particleidx]
  # copy!(fr.Y, fr.X[:,fr.gwp.particleidx])
  for i in 1:fr.zDim
    fr.Y[fr.p[i]] = X[i]
  end
  nothing
end

function (p::GenericWrapParam)(res::Vector{Float64}, x::Vector{Float64})

  # TODO -- move to inner lambda that is defined once against p.params...
  # approximates by not considering cross indices among parameters
  # @show length(p.params), p.varidx, p.particleidx, size(x), size(res), size(p.measurement)
  p.params[p.varidx][:, p.particleidx] = x
  # p.usrfnc!(res, p.particleidx, p.measurement, p.params...)
  # who are active hypotheses?  p.params[p.activehypo]...
  p.usrfnc!(res, p.factormetadata, p.particleidx, p.measurement, view(p.params,p.activehypo)...)
  # p.usrfnc!(res, p.factormetadata, p.particleidx, p.measurement, p.params[p.activehypo]...)
end

function (fr::FastRootGenericWrapParam)( res::Vector{Float64}, x::Vector{Float64} )


  shuffleXAltD!(fr, x)
  fr.gwp( res, fr.Y )
end

function numericRootGenericRandomizedFnc!(
            frl::FastRootGenericWrapParam{T};
            perturb::Float64=1e-10,
            testshuffle::Bool=false ) where {T <: FunctorPairwiseMinimize}
  #

  fill!(frl.res, 0.0) # 1:frl.xDim
  r = optimize( frl.gg, frl.X[:, frl.gwp.particleidx] ) # frl.gg
  # TODO -- clearly lots of optmization to be done here
  frl.Y[:] = r.minimizer
  frl.X[:,frl.gwp.particleidx] = frl.Y
  nothing
end

function numericRootGenericRandomizedFnc!(
            fr::FastRootGenericWrapParam{T};
            perturb::Float64=1e-10,
            testshuffle::Bool=false ) where {T <: FunctorPairwise}
  #
  # info("numericRootGenericRandomizedFnc! FastRootGenericWrapParam{$T}")
  # @show fr.zDim, fr.xDim, fr.gwp.partial, fr.gwp.particleidx

  if fr.zDim < fr.xDim && !fr.gwp.partial || testshuffle
    # less measurement dimensions than variable dimensions -- i.e. shuffle
    shuffle!(fr.p)
    for i in 1:fr.xDim
      fr.perturb[1:fr.zDim] = perturb*randn(fr.zDim)
      fr.X[fr.p[1:fr.zDim], fr.gwp.particleidx] += fr.perturb
      r = nlsolve(  fr,
                    fr.X[fr.p[1:fr.zDim], fr.gwp.particleidx] # this is x0
                 )
      if r.f_converged
        shuffleXAltD!( fr, r.zero )
        break;
      else
        # TODO -- report on this bottleneck, useful for optimization of code
        # @show i, fr.p, fr.xDim, fr.zDim
        temp = fr.p[end]
        fr.p[2:end] = fr.p[1:(end-1)]
        fr.p[1] = temp
        if i == fr.xDim
          error("numericRootGenericRandomizedFnc could not converge, i=$(i), fr.gwp.usrfnc!=$(typeof(fr.gwp.usrfnc!))")
        end
      end
    end
    #shuffleXAltD!( fr, r.zero ) # moved up
  elseif fr.zDim >= fr.xDim && !fr.gwp.partial
    # equal or more measurement dimensions than variable dimensions -- i.e. don't shuffle
    fr.perturb[1:fr.xDim] = perturb*randn(fr.xDim)
    fr.X[1:fr.xDim, fr.gwp.particleidx] += fr.perturb[1:fr.xDim] # moved up
    r = nlsolve( fr.gwp, fr.X[1:fr.xDim,fr.gwp.particleidx] )
    if sum(isnan.(( r ).zero)) == 0
      fr.Y[1:fr.xDim] = ( r ).zero
    else
      warn("got NaN, fr.gwp.particleidx = $(fr.gwp.particleidx), r=$(r)")
      @show fr.gwp.usrfnc!
      for thatlen in 1:length(fr.gwp.params)
        @show thatlen, fr.gwp.params[thatlen][:, fr.gwp.particleidx]
      end
    end
  elseif fr.gwp.partial
    # improve memory management in this function
    fr.p[1:length(fr.gwp.usrfnc!.partial)] = Int[fr.gwp.usrfnc!.partial...] # TODO -- move this line up and out of inner loop
    r = nlsolve(  fr,
                  fr.X[fr.p[1:fr.zDim], fr.gwp.particleidx] # this is x0
               )
    shuffleXAltD!( fr, r.zero )
  else
    error("Unresolved numeric solve case")
  end
  fr.X[:,fr.gwp.particleidx] = fr.Y
  nothing
end



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
