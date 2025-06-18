# TODO: KEEP
# function compareAll(Al::T1, Bl::T2; show::Bool=true, skip::Vector{Symbol}=Symbol[])::Bool where {T1 <: Union{SingleThreaded, MultiThreaded}, T2 <: Union{SingleThreaded, MultiThreaded}}
#   return T1 == T2
# end
import DistributedFactorGraphs: compare

# These are now moved to DistributedFactorGraphs, with the exceptions of
# the functions with IIF-specific parameters.
# To extend these, import the relevant DFG compareX function and overload it.

function Base.isapprox(
  p1::Union{<:BallTreeDensity, <:ManifoldKernelDensity},
  p2::Union{<:BallTreeDensity, <:ManifoldKernelDensity};
  atol::Real = 1e-6,
)
  #
  return mmd(p1, p2) < atol
end



## FIXME, FIGURE OUT HOW TO DEPRECATE BELOW ==============================================

function compareAllSpecial(
  A::T1,
  B::T2;
  skip = Symbol[],
  show::Bool = true,
) where {T1 <: CommonConvWrapper, T2 <: CommonConvWrapper}
  #
  if T1 != T2
    @warn "CCW types T1 and T2 not equal=>" T1 T2
    # return false
  end
  # FIXME still issues with compare, skipping :vartypes https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/434
  return compareAll(A, B; skip = union(skip, [:vartypes]), show = show)
end

function compare(
  p1::Union{<:BallTreeDensity, <:ManifoldKernelDensity},
  p2::Union{<:BallTreeDensity, <:ManifoldKernelDensity},
)
  #
  return compareAll(p1.bt, p2.bt; skip = [:calcStatsHandle; :data]) &&
         compareAll(p1, p2; skip = [:calcStatsHandle; :bt])
end

function compare(c1::TreeClique, c2::TreeClique)
  #
  TP = true
  TP = TP && c1.id == c2.id
  # data
  @warn "skipping ::TreeClique compare of data"
  # TP = TP && compare(c1.data, c2.data)

  # attributes
  @warn "only comparing keys of TreeClique attributes"
  TP = TP && collect(keys(c1.attributes)) == collect(keys(c2.attributes))

  return TP
end

#
