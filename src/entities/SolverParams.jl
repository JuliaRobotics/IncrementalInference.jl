

"""
$(TYPEDEF)

Solver parameters for the DistributedFactoGraph.

Dev Notes
- FIXME change to using kwargs from Parameters.jl
- TODO remove NothingUnion
- TODO Upgrade to common @kwargs struct approach
"""
mutable struct SolverParams <: DFG.AbstractParams
  dimID::Int
  # TODO remove NothingUnion
  registeredModuleFunctions::NothingUnion{Dict{Symbol, Function}} # remove from
  reference::NothingUnion{Dict{Symbol, Tuple{Symbol, Vector{Float64}}}}
  stateless::Bool
  qfl::Int # Quasi fixed length
  isfixedlag::Bool # true when adhering to qfl window size for solves
  limitfixeddown::Bool # if true, then fixed lag will also not update marginalized during down (default false)
  # new functions
  incremental::Bool
  useMsgLikelihoods::Bool
  upsolve::Bool
  downsolve::Bool
  drawtree::Bool
  drawCSMIters::Bool
  showtree::Bool
  drawtreerate::Float64
  dbg::Bool
  async::Bool
  limititers::Int
  N::Int
  multiproc::Bool
  logpath::String
  graphinit::Bool
  treeinit::Bool # still experimental with known errors
  limittreeinit_iters::Int
  algorithms::Vector{Symbol} # list of algorithms to run [:default] is mmisam
  spreadNH::Float64 # experimental, entropy spread adjustment used for both null hypo cases.
  inflation::Float64 # experimental, how much to disperse particles before convolution solves, #1051
  inflateCycles::Int
  gibbsIters::Int
  maxincidence::Int # maximum incidence to a variable in an effort to enhance sparsity
  alwaysFreshMeasurements::Bool
  # should factor gradients be calculated or attempted (UNDER DEVELOPMENT, 21Q3)
  attemptGradients::Bool
  devParams::Dict{Symbol,String}
  #
end

SolverParams(;dimID::Int=0,
              registeredModuleFunctions=nothing,
              reference=nothing,
              stateless::Bool=false,
              qfl::Int=99999999999,
              isfixedlag::Bool=false,
              limitfixeddown::Bool=false,
              incremental::Bool=true,
              useMsgLikelihoods::Bool=false,
              upsolve::Bool=true,
              downsolve::Bool=true,
              drawtree::Bool=false,
              drawCSMIters::Bool=true,
              showtree::Bool=false,
              drawtreerate::Float64=0.5,
              dbg::Bool=false,
              async::Bool=false,
              limititers::Int=500,
              N::Int=100,
              multiproc::Bool=1 < nprocs(),
              logpath::String="/tmp/caesar/$(now())",
              graphinit::Bool=true,
              treeinit::Bool=false,
              limittreeinit_iters::Int=10,
              algorithms::Vector{Symbol}=[:default],
              spreadNH::Real=3.0,
              inflation::Real=5.0,
              inflateCycles::Int=3,
              gibbsIters::Int=3,
              maxincidence::Int=500,
              alwaysFreshMeasurements::Bool=true,
              attemptGradients::Bool=true,
              devParams::Dict{Symbol,String}=Dict{Symbol,String}()
            ) = begin useMsgLikelihoods==true && @warn "useMsgLikelihoods is under development, use with care, see #1010"
                SolverParams( dimID,
                              registeredModuleFunctions,
                              reference,
                              stateless,
                              qfl,
                              isfixedlag,
                              limitfixeddown,
                              incremental,
                              useMsgLikelihoods,
                              upsolve,
                              downsolve,
                              drawtree,
                              drawCSMIters,
                              showtree,
                              drawtreerate,
                              dbg,
                              async,
                              limititers,
                              N,
                              multiproc,
                              logpath,
                              graphinit,
                              treeinit,
                              limittreeinit_iters,
                              algorithms,
                              spreadNH,
                              inflation,
                              inflateCycles,
                              gibbsIters,
                              maxincidence,
                              alwaysFreshMeasurements,
                              attemptGradients,
                              devParams )
            end
#

