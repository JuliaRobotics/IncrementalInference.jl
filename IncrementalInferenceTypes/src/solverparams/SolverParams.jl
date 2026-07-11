
"""
$(TYPEDEF)

Solver parameters for the DistributedFactoGraph.

Dev Notes
- FIXME change to using kwargs from Parameters.jl
- TODO remove NothingUnion
- TODO Upgrade to common @kwargs struct approach
"""
Base.@kwdef mutable struct SolverParams <: DFG.AbstractDFGParams
  dimID::Int = 0
  reference::Union{Nothing, Dict{Symbol, Tuple{Symbol, Vector{Float64}}}} = nothing
  stateless::Bool = false
  """ Quasi fixed length """
  qfl::Int = (2^(Sys.WORD_SIZE - 1) - 1)
  """ true when adhering to qfl window size for solves """
  isfixedlag::Bool = false          
  """ if true, then fixed lag will not update marginalized during down pass on tree """
  limitfixeddown::Bool = false      
  """ use incremental tree updates, TODO consolidate with recycling """
  incremental::Bool = true          
  """ Experimental, insert differential factors from upward joints """
  useMsgLikelihoods::Bool = false   
  """ do tree upsolve """
  upsolve::Bool = true              
  """ do tree downsolve """
  downsolve::Bool = true            
  """ draw tree during solve """
  drawtree::Bool = false            
  """ show CSM iteration count on tree visualization """
  drawCSMIters::Bool = true         
  showtree::Bool = false
  """ how fast should the tree vis file be redrawn """
  drawtreerate::Float64 = 0.5       
  """ Experimental, enable additional tier debug features """
  dbg::Bool = false                 
  """ do not block on CSM tasks """
  async::Bool = false               
  """ limit number of steps CSMs can take """
  limititers::Int = 500             
  """ default number of particles """
  N::Int = 100                      
  """ should Distributed.jl tree solve compute features be used """
  multiproc::Bool = 1 < nprocs()    
  """ "/tmp/caesar/logs/$(now())" # unique temporary file storage location for a solve """
  logpath::String = joinpath(tempdir(),"caesar","logs","$(now(DFG.UTC))") 
  """ default to graph-based initialization of variables """
  graphinit::Bool = true            
  """ init variables on the tree """
  treeinit::Bool = false             
  limittreeinit_iters::Int = 10
  """ list of algorithms to run [:default] is mmisam """
  algorithms::Vector{Symbol} = [:default, :parametric] 
  """ entropy spread adjustment used for both null hypo cases. """
  spreadNH::Float64 = 3.0           
  """ how much to disperse particles before convolution solves, #1051 """
  inflation::Float64 = 5.0          
  """ minimum nullhypo for relative factors sibling to multihypo factors onto a specific variable. """
  nullSurplusAdd::Float64 = 0.3     
  """ repeat convolutions for inflation to occur """
  inflateCycles::Int = 3            
  """ number of Gibbs cycles to take per clique iteration variables """
  gibbsIters::Int = 3               
  """ maximum incidence to a variable in an effort to enhance sparsity """
  maxincidence::Int = 500           
  """ Development feature on whether new samples should be sampled at each Gibbs cycle convolution """
  alwaysFreshMeasurements::Bool = true 
  """ should factor gradients be calculated or attempted (UNDER DEVELOPMENT, 21Q3) """
  attemptGradients::Bool = false    
  """ empty container for new features, allowing workaround for breaking changes and legacy """
  devParams::Dict{Symbol, String} = Dict{Symbol, String}() 
  #
end

StructTypes.omitempties(::Type{SolverParams}) = (:reference,)

#
Base.convert(::Type{SolverParams}, ::DFG.NoSolverParams) = begin 
  @warn "FIXME Why converting NoSolverParams to SolverParams?"
  SolverParams()
end

