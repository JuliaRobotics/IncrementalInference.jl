

"""
$(TYPEDEF)

Solver parameters for the DistributedFactoGraph.

Dev Notes
- FIXME change to using kwargs from Parameters.jl
- TODO remove NothingUnion
- TODO Upgrade to common @kwargs struct approach
"""
Base.@kwdef mutable struct SolverParams <: DFG.AbstractParams
  dimID::Int = 0
  reference::NothingUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}} = nothing
  stateless::Bool = false
  qfl::Int = 99999999999            # Quasi fixed length
  isfixedlag::Bool = false          # true when adhering to qfl window size for solves
  limitfixeddown::Bool = false      # if true, then fixed lag will also not update marginalized during down (default false)
  incremental::Bool = true          # use incremental tree updates, TODO consolidate with recycling
  useMsgLikelihoods::Bool = false   # Experimental, insert differential factors from upward joints
  upsolve::Bool = true              # do tree upsolve
  downsolve::Bool = true            # do tree downsolve
  drawtree::Bool = false            # draw tree during solve
  drawCSMIters::Bool = true         # show CSM iteration count on tree visualization
  showtree::Bool = false
  drawtreerate::Float64 = 0.5       # how fast should the tree vis file be redrawn
  dbg::Bool = false                 # Experimental, enable additional tier debug features
  async::Bool = false               # do not block on CSM tasks
  limititers::Int = 500             # limit number of steps CSMs can take
  N::Int = 100                      # default number of particles
  multiproc::Bool = 1 < nprocs()    # should Distributed.jl tree solve compute features be used
  logpath::String = "/tmp/caesar/$(now())" # unique temporary file storage location for a solve
  graphinit::Bool = true            # default to graph-based initialization of variables
  treeinit::Bool =false             # Experimental, init variables on the tree
  limittreeinit_iters::Int = 10
  algorithms::Vector{Symbol} = [:default] # list of algorithms to run [:default] is mmisam
  spreadNH::Float64 = 3.0           # Experimental, entropy spread adjustment used for both null hypo cases.
  inflation::Float64 = 5.0          # Experimental, how much to disperse particles before convolution solves, #1051
  inflateCycles::Int = 3            # repeat convolutions for inflation to occur
  gibbsIters::Int = 3               # number of Gibbs cycles to take per clique iteration variables
  maxincidence::Int = 500           # maximum incidence to a variable in an effort to enhance sparsity
  alwaysFreshMeasurements::Bool = true # Development feature on whether new samples should be sampled at each Gibbs cycle convolution
  attemptGradients::Bool = false    # should factor gradients be calculated or attempted (UNDER DEVELOPMENT, 21Q3)
  devParams::Dict{Symbol,String} = Dict{Symbol,String}() # empty container for new features, allowing workaround for breaking changes and legacy
  #
end



#