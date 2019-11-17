# TODO: KEEP
# function compareAll(Al::T1, Bl::T2; show::Bool=true, skip::Vector{Symbol}=Symbol[])::Bool where {T1 <: Union{SingleThreaded, MultiThreaded}, T2 <: Union{SingleThreaded, MultiThreaded}}
#   return T1 == T2
# end
import DistributedFactorGraphs: compare, compareAllSpecial

# These are now moved to DistributedFactorGraphs, with the exceptions of
# the functions with IIF-specific parameters.
# To extend these, import the relevant DFG compareX function and overload it.

function compare(p1::BallTreeDensity, p2::BallTreeDensity)::Bool
  return compareAll(p1.bt,p2.bt, skip=[:calcStatsHandle; :data]) &&
         compareAll(p1,p2, skip=[:calcStatsHandle; :bt])
end

function compareAllSpecial(A::T1, B::T2;
                    skip=Symbol[], show::Bool=true) where {T1 <: CommonConvWrapper, T2 <: CommonConvWrapper}
  #
  T1 != T2 && return false
  return compareAll(A, B, skip=skip, show=show)
end
