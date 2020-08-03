


##==============================================================================
## Delete at end v0.14.x
##==============================================================================


export Pairwise, Singleton
export FunctorSingletonNH, FunctorPairwiseNH
export FunctorInferenceType, FunctorPairwise, FunctorPairwiseMinimize, FunctorSingleton



# TODO been replaced by Functor types, but may be reused for non-numerical cases
abstract type Pairwise <: InferenceType end
abstract type Singleton <: InferenceType end

# TODO deprecate with standard null hypothesis only
abstract type FunctorSingletonNH <: FunctorSingleton end
abstract type FunctorPairwiseNH <: FunctorPairwise end
# abstract type FunctorPairwiseNHMinimize <: FunctorPairwiseMinimize end # TODO


# const FGG = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Array{Graphs.ExVertex,1},Array{Array{Graphs.Edge{Graphs.ExVertex},1},1}}
# const FGGdict = Graphs.GenericIncidenceList{Graphs.ExVertex,Graphs.Edge{Graphs.ExVertex},Dict{Int,Graphs.ExVertex},Dict{Int,Array{Graphs.Edge{Graphs.ExVertex},1}}}


"""
    $(SIGNATURES)

Do true and null hypothesis computations based on data structures prepared earlier -- specific to `FunctorPairwiseNH`.  This function will be merged into a standard case for `AbstractRelativeFactor/Minimize` in the future.
"""
function computeAcrossNullHypothesis!(ccwl::CommonConvWrapper{T},
                                      allelements,
                                      nhc,
                                      ENT  ) where {T <: FunctorPairwiseNH}
  #
  # TODO --  Threads.@threads see area4 branch
  for n in allelements
    # ccwl.gwp(x, res)
    if nhc[n] != 0
      ccwl.cpt[Threads.threadid()].particleidx = n
      numericRootGenericRandomizedFnc!( ccwl )
    else
      ccwl.params[ccwl.varidx][:,n] += rand(ENT)
    end
  end
  nothing
end


"""
    $(SIGNATURES)

Prepare data required for null hypothesis cases during convolution.
"""
function assembleNullHypothesis(ccwl::CommonConvWrapper{T},
                                maxlen::Int,
                                spreadfactor::Real=10 ) where {T}
  #
  @warn "this assembleNullHypothesis method has been deprecated for e.g. `addFactor!(; nullhypo=0.1)` instead."
  nhc = rand(ccwl.usrfnc!.nullhypothesis, maxlen) .- 1
  arr = ccwl.params[ccwl.varidx]
  ENT = generateNullhypoEntropy(arr, maxlen, spreadfactor)
  allelements = 1:maxlen
  return allelements, nhc, ENT
end


