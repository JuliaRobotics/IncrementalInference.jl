


"""
$(TYPEDEF)

Main factor memory container used during inference operations -- i.e. values specific to one complete convolution operation

Notes
- CCW does not get serialized / persisted
- At writing, the assumption is there is just one CCW per factor
- Any multithreaded design needs to happens as sub-constainers inside CCW or otherwise, to carry separate memory.
- Since #467, `CalcFactor` is the only type 'seen by user' during `getSample` or function residual calculations `(cf::CalcFactor{<:MyFactor})`, s.t. `MyFactor <: AbstractRelative___`
- There also exists a `CalcFactorMahalanobis` for parameteric computations using as much of the same mechanics as possible.
- CCW is consolidated object of other previous types, FMD CPT CF CFM.

Related 

[`CalcFactor`](@ref), [`CalcFactorMahalanobis`](@ref)
"""
Base.@kwdef struct CommonConvWrapper{
  T <: AbstractFactor, 
  VT <: Tuple,
  TP <: Base.RefValue{<:Tuple},
  CT,
  AM <: AbstractManifold,
  HR <: HypoRecipeCompute,
  MT, 
  G,
  KCF <: Union{Nothing, <: Channel{<: CalcFactor}}
} <: FactorCache
  # Basic factor topological info
  """ Values consistent across all threads during approx convolution """
  usrfnc!::T # user factor / function
  """ Ordered tuple of all variables connected to this factor """
  fullvariables::VT
  # shortcuts to numerical containers
  """ Numerical containers for all connected variables.  Hypo selection needs to be passed 
      to each hypothesis evaluation event on user function via CalcFactor, #1321.
      Points directly at the variable VND.val (not a deepcopy). """
  varValsAll::TP
  """ dummy cache value to be deep copied later for each of the CalcFactor instances """
  dummyCache::CT = nothing
  # derived config parameters for this factor
  """ Factor manifold definition for frequent use (not the variables manifolds) """
  manifold::AM = getManifold(usrfnc!)
  """ Which dimensions does this factor influence.  Sensitive (mutable) to both which 'solvefor index' variable and whether the factor is partial dimension """
  partialDims::Vector{<:Integer} = collect(1:manifold_dimension(manifold))
  """ is this a partial constraint as defined by the existance of factor field `.partial::Tuple` """
  partial::Bool = false
  """ probability that this factor is wholly incorrect and should be ignored during solving """
  nullhypo::Float64 = 0.0
  """ inflationSpread particular to this factor (by how much to dispurse the belief initial values before numerical optimization is run).  Analogous to stochastic search """
  inflation::Float64 = SolverParams().inflation
  """ multihypo specific field containers for recipe of hypotheses to compute """
  hyporecipe::HR = HypoRecipeCompute(;activehypo=collect(1:length(varValsAll)))
  # buffers and indices to point numerical computations to specific memory locations
  """ user defined measurement values for each approxConv operation
      FIXME make type stable, JT should now be type stable if rest works.
      SUPER IMPORTANT, if prior=>point or relative=>tangent, see #1661 
      can be a Vector{<:Tuple} or more direct Vector{<: pointortangenttype} """
  measurement::Vector{MT} = Vector(Vector{Float64}())
  """ which index is being solved for in params? """
  varidx::Base.RefValue{Int} = Ref(1)
  """ Consolidation from CPT, the actual particle being solved at this moment """
  particleidx::Base.RefValue{Int} = Ref(1)
  """ working memory to store residual for optimization routines """
  res::Vector{Float64} = zeros(manifold_dimension(manifold))
  """ experimental feature to embed gradient calcs with ccw """
  _gradients::G = nothing
    """ working memory to store residual from optimization routines """
  keepCalcFactor::KCF
end


#
