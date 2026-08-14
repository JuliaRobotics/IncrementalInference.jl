## ======================================================================================
## Density bundles: the types, and the layout accessors that come with them.
##
## Nothing here touches a graph and nothing *does* anything to a belief.  The algebra is in
## `parametric/services/BeliefAlgebra.jl`; the graph boundary in
## `parametric/services/ParametricTreeSolve.jl`.
##
## ======================================================================================

"""
    $TYPEDEF

A Gaussian in **canonical** (information) form over stacked tangent coordinates `ξ`, with precision `Λ`
and information vector `η = Λμ`:

    φ(ξ) = ½ ξᵀ Λ ξ - ηᵀ ξ - logmass,      p(ξ) ∝ exp(-φ(ξ))

*Co*tangent because that is where these parameters live — `η` is a covector, `Λ` a bilinear form on the
tangent space.  [`TangentNormal`](@ref) is the moment-form dual.

`η === nothing` is the **centred** case `μ = 0` ([`isCentred`](@ref)); `η` present is *extended*.

$(TYPEDFIELDS)
"""
struct CotangentNormal{ET <: Union{Nothing, Vector{Float64}}}
  """precision (information matrix); may be legitimately singular"""
  Λ::Matrix{Float64}
  """information vector `Λμ`, or `nothing` for the centred case `μ = 0`"""
  η::ET
  """log of the mass, entering `φ` **negated** so that `p ∝ exp(+logmass)` — more `logmass` is more
  mass.  Gauss-Newton stores `−½‖r‖²`, a worse fit being less mass."""
  logmass::Float64

  function CotangentNormal(
    Λ::AbstractMatrix,
    η::Union{Nothing, AbstractVector} = nothing,
    logmass::Real = 0.0,
  )
    size(Λ, 1) == size(Λ, 2) ||
      error("CotangentNormal: Λ is $(size(Λ,1))×$(size(Λ,2)), must be square")
    isnothing(η) || length(η) == size(Λ, 1) ||
      error("CotangentNormal: η has length $(length(η)) but Λ spans $(size(Λ,1)) dimensions")
    η_ = isnothing(η) ? nothing : Vector{Float64}(η)
    return new{typeof(η_)}(Matrix{Float64}(Λ), η_, Float64(logmass))
  end
end

"""
    $SIGNATURES

Is the mean zero in tangent coordinates.
"""
isCentred(fibre::CotangentNormal) = isnothing(fibre.η)

"""
    $SIGNATURES
Coordinate dimension the density spans.
"""
Base.length(fibre::CotangentNormal) = size(fibre.Λ, 1)

"""
    $SIGNATURES
Information vector, materializing zeros for the zero mean case so callers need not branch.
"""
getInfovector(fibre::CotangentNormal) = @something(fibre.η, zeros(length(fibre)))


"""
    $TYPEDEF

A Gaussian in **moment** form over stacked tangent coordinates — the dual of [`CotangentNormal`](@ref),
with mean `μ` (a tangent vector) and covariance `Σ`.

$(TYPEDFIELDS)
"""
struct TangentNormal{MT <: Union{Nothing, Vector{Float64}}}
  """covariance; on the observable subspace only, where the fibre came from a singular precision"""
  Σ::Matrix{Float64}
  """mean tangent vector, or `nothing` for the centred case `μ = 0`"""
  μ::MT
  """log of the unnormalized integral — see [`CotangentNormal`](@ref)"""
  logmass::Float64

  function TangentNormal(
    Σ::AbstractMatrix,
    μ::Union{Nothing, AbstractVector} = nothing,
    logmass::Real = 0.0,
  )
    size(Σ, 1) == size(Σ, 2) || error("TangentNormal: Σ is $(size(Σ,1))×$(size(Σ,2)), must be square")
    isnothing(μ) || length(μ) == size(Σ, 1) ||
      error("TangentNormal: μ has length $(length(μ)) but Σ spans $(size(Σ,1)) dimensions")
    μ_ = isnothing(μ) ? nothing : Vector{Float64}(μ)
    return new{typeof(μ_)}(Matrix{Float64}(Σ), μ_, Float64(logmass))
  end
end

isCentred(fibre::TangentNormal) = isnothing(fibre.μ)
Base.length(fibre::TangentNormal) = size(fibre.Σ, 1)

# ======================================================================================
# The joint layout, and points of a density bundle over it.
# ======================================================================================

"""
    $TYPEDEF

One block of the base manifold: `count` consecutive variables that all share state type `statetype`.

$(TYPEDFIELDS)
"""
struct PowerStateType
  """state type shared by every variable in the block"""
  statetype::DataType
  """how many variables the block holds"""
  count::Int
end

Base.length(b::PowerStateType) = b.count
getDimension(b::PowerStateType) = getDimension(b.statetype)

getManifold(b::PowerStateType) = NPowerManifold(getManifold(b.statetype), b.count)

"""
    $SIGNATURES
Total stacked coordinate dimension of a block.
"""
getCoordwidth(b::PowerStateType) = b.count * getDimension(b.statetype)

"""
    $TYPEDEF

Which variables a joint density spans, what geometry each has, and which stacked coordinates belong to
which.

Separated out because it is the part that does *not* vary: [`pullback`](@ref), [`fuse`](@ref),
[`stepBelief`](@ref) and [`calcMomentform`](@ref) all hand the *same* layout object on, and every
component of a mixture over one variable set will share it too.  Rebuilding it would allocate an
identical `ranges` and `OrderedDict` per clique, per message, per sweep.

`labels` is grouped by state type — the order `getVarIntLabelMap` produces — so blocks are contiguous
and any order-preserving subset is still grouped, which is what lets [`subsetLayout`](@ref) work
without a graph.  The geometry lives in `blocks` ([`PowerStateType`](@ref)), one entry per block rather
than one per variable; `blockof` maps back the other way for the per-variable views.

$(TYPEDFIELDS)
"""
struct JointLayout
  """variables whose tangent coordinates are stacked, in order"""
  labels::Vector{Symbol}
  """the base manifold's power blocks, in coordinate order — the geometry, and the block structure the
  base point must match ([`DensityBundlePoint`](@ref) checks it)"""
  blocks::Vector{PowerStateType}
  """which block each entry of `labels` falls in; the inverse of `blocks`' runs"""
  blockof::Vector{Int}
  """index of the FIRST entry of `labels` in each block — so a block-wise walk needs no running offset"""
  blockstarts::Vector{Int}
  """stacked coordinate range of each entry of `labels`"""
  ranges::Vector{UnitRange{Int}}
  """`labels[i] => i`, so a receiver can address a sender's variables without sharing its layout."""
  index::OrderedDict{Symbol, Int}

  function JointLayout(
    labels::AbstractVector{Symbol},
    statetypes::AbstractVector,
    blocklengths::Union{Nothing, AbstractVector{Int}} = nothing,
  )
    length(statetypes) == length(labels) ||
      error("JointLayout: $(length(statetypes)) state types for $(length(labels)) labels")
    ranges = UnitRange{Int}[]
    start = 1
    for T in statetypes
      d = getDimension(T)
      push!(ranges, start:(start + d - 1))
      start += d
    end
    # Blocking is GIVEN when the caller has a partition the types cannot express — a frontal/separator
    # split puts the same state type on both sides, and consecutive runs would silently merge the two
    # into one block that no longer matches the base point's container.  Otherwise it is inferred as
    # runs of consecutive equal state types, the grouping `getVarIntLabelMap` produces.
    blocks = PowerStateType[]
    blockof = Vector{Int}(undef, length(statetypes))
    blockstarts = Int[]
    if isnothing(blocklengths)
      for (k, T) in enumerate(statetypes)
        if k == 1 || T !== statetypes[k - 1]
          push!(blocks, PowerStateType(T, 1))
          push!(blockstarts, k)
        else
          blocks[end] = PowerStateType(T, blocks[end].count + 1)
        end
        blockof[k] = length(blocks)
      end
    else
      sum(blocklengths) == length(labels) || error(
        "JointLayout: blocklengths sum to $(sum(blocklengths)) but there are $(length(labels)) labels",
      )
      k = 1
      for len in blocklengths
        len > 0 || error("JointLayout: empty block")
        T = statetypes[k]
        push!(blocks, PowerStateType(T, len))
        push!(blockstarts, k)
        for j in k:(k + len - 1)
          statetypes[j] === T ||
            error("JointLayout: block starting at $(labels[k]) mixes $(T) and $(statetypes[j])")
          blockof[j] = length(blocks)
        end
        k += len
      end
    end
    return new(
      collect(Symbol, labels),
      blocks,
      blockof,
      blockstarts,
      ranges,
      OrderedDict{Symbol, Int}(zip(labels, eachindex(labels))),
    )
  end
end

"""
    $SIGNATURES
State type of the `k`-th variable — the per-variable view of [`PowerStateType`](@ref) `blocks`.
"""
getStatetype(layout::JointLayout, k::Int) = layout.blocks[layout.blockof[k]].statetype

"""
    $SIGNATURES
How many variables fall in each block.
"""
getBlocklengths(layout::JointLayout) = Int[b.count for b in layout.blocks]

"""
    $SIGNATURES
Flat per-variable state types — the denormalized view of [`PowerStateType`](@ref) `blocks`.

"""
getStatetypes(layout::JointLayout) =
  DataType[getStatetype(layout, k) for k in eachindex(layout.labels)]

"""
    $TYPEDEF

A point of a **density bundle**: base points on a (product) manifold together with the density in the
fibre above them, over a shared [`JointLayout`](@ref).  With a [`CotangentNormal`](@ref) fibre this is
the exponentially wrapped Gaussian `q = exp(p, ε)`, `ε ~ N(μ, Λ⁻¹)`.

The bundle is **trivial** — `G × (ℝ × ℝⁿ × Sym⁺(n))` for a Lie group `G` — which is why the fibre can be
a type of its own with no reference to the base.  The *transport* between fibres is not, see [`pullback`](@ref).

`point`: an `ArrayPartition` of per-variable points is one point of the product manifold.
Slice it against the joint precision by label — `getBasepoint(bundle, s)` and
`getCoordrange(bundle, s)` — or walk both with [`eachvariable`](@ref).

$(TYPEDFIELDS)
"""
struct DensityBundlePoint{PT, FT}
  """which variables, what geometry, which coordinates — shared, see [`JointLayout`](@ref)"""
  layout::JointLayout
  """the base point on the product manifold, blocked by state type (an `ArrayPartition`, as the
  parametric solvers use)"""
  point::PT
  """the density over the stacked tangent coordinates at `point`"""
  fibre::FT

  function DensityBundlePoint(layout::JointLayout, point, fibre)
    # The block structure is a FUNCTION of the layout, not free — this check is what lets
    # `getBasepoint` index the `ArrayPartition` flatly.  A total-count check would not do: a container
    # blocked differently would pass it, leaving `getBasepoint` and `getCoordrange` disagreeing about
    # which variable is which.
    length(point.x) == length(layout.blocks) || error(
      "DensityBundlePoint: base point has $(length(point.x)) blocks but $(layout.labels) group into $(length(layout.blocks))",
    )
    for (b, block) in enumerate(point.x)
      length(block) == layout.blocks[b].count || error(
        "DensityBundlePoint: block $b holds $(length(block)) points, layout expects $(layout.blocks[b].count)",
      )
    end
    length(fibre) == length(layout) || error(
      "DensityBundlePoint: fibre spans $(length(fibre)) but $(layout.labels) span $(length(layout)) dimensions",
    )
    return new{typeof(point), typeof(fibre)}(layout, point, fibre)
  end
end

"""
    $SIGNATURES
Build the layout in place, for a bundle point that does not share one with anything.
"""
DensityBundlePoint(
  point,
  labels::AbstractVector{Symbol},
  statetypes::AbstractVector,
  fibre,
) = DensityBundlePoint(JointLayout(labels, statetypes), point, fibre)

"""
    $SIGNATURES
Build the [`CotangentNormal`](@ref) fibre in place, for the common canonical-form case.
"""
DensityBundlePoint(
  layout::JointLayout,
  point,
  Λ::AbstractMatrix,
  η::Union{Nothing, AbstractVector} = nothing,
  logmass::Real = 0.0,
) = DensityBundlePoint(layout, point, CotangentNormal(Λ, η, logmass))

DensityBundlePoint(
  point,
  labels::AbstractVector{Symbol},
  statetypes::AbstractVector,
  Λ::AbstractMatrix,
  η::Union{Nothing, AbstractVector} = nothing,
  logmass::Real = 0.0,
) = DensityBundlePoint(JointLayout(labels, statetypes), point, CotangentNormal(Λ, η, logmass))

# Layout questions are answered by the bundle point directly: storage lives in `layout`, the interface
# does not move.
@inline function Base.getproperty(bundle::DensityBundlePoint, s::Symbol)
  s === :labels && return getfield(bundle, :layout).labels
  # derived, not stored — the geometry lives per BLOCK now, see [`PowerStateType`](@ref).  Kept for
  # the flat per-variable views; a hot walk should go through `layout.blocks` instead.
  s === :statetypes && return getStatetypes(getfield(bundle, :layout))
  s === :blocks && return getfield(bundle, :layout).blocks
  s === :ranges && return getfield(bundle, :layout).ranges
  s === :index && return getfield(bundle, :layout).index
  return getfield(bundle, s)
end

Base.propertynames(::DensityBundlePoint) =
  (:layout, :point, :fibre, :labels, :statetypes, :blocks, :ranges, :index)

"""
    $SIGNATURES
Total stacked coordinate dimension — a property of the layout, so it holds whatever the fibre is.
"""
Base.length(layout::JointLayout) = isempty(layout.ranges) ? 0 : last(layout.ranges[end])
Base.length(bundle::DensityBundlePoint) = length(bundle.layout)

isCentred(bundle::DensityBundlePoint) = isCentred(bundle.fibre)
getInfovector(bundle::DensityBundlePoint) = getInfovector(bundle.fibre)

"""
    $SIGNATURES
Does this cover variable `label`?.
"""
hasLabel(layout::JointLayout, label::Symbol) = haskey(layout.index, label)
hasLabel(bundle::DensityBundlePoint, label::Symbol) = hasLabel(bundle.layout, label)

"""
    $SIGNATURES
Base point of variable `label`.
"""
getBasepoint(bundle::DensityBundlePoint, label::Symbol) = bundle.point[bundle.index[label]]

"""
    $SIGNATURES
Stacked coordinate range of variable `label` — slices `Λ`, `η`, and any covariance built from them.
"""
getCoordrange(layout::JointLayout, label::Symbol) = layout.ranges[layout.index[label]]
getCoordrange(bundle::DensityBundlePoint, label::Symbol) = getCoordrange(bundle.layout, label)

"""
    $SIGNATURES
Tangent dimension of each entry of `labels`.
"""
getDimensions(layout::JointLayout) = Int[length(r) for r in layout.ranges]
getDimensions(bundle::DensityBundlePoint) = getDimensions(bundle.layout)

"""
    $SIGNATURES
Stacked coordinate indices of `labels`, in the order given.
"""
function getCoordindices(layout::JointLayout, labels::AbstractVector{Symbol})
  for label in labels
    hasLabel(layout, label) || error("JointLayout has no coordinates for $label")
  end
  return reduce(vcat, (collect(getCoordrange(layout, label)) for label in labels); init = Int[])
end

getCoordindices(bundle::DensityBundlePoint, labels::AbstractVector{Symbol}) =
  getCoordindices(bundle.layout, labels)

"""
    $SIGNATURES

Coordinate indices of `labels` as a **`UnitRange`** when they are contiguous and in layout order, and
as the explicit index vector otherwise.
"""
function getCoordspan(layout::JointLayout, labels::AbstractVector{Symbol})
  isempty(labels) && return 1:0
  # Decided on the LABEL indices, so the contiguous case never materializes an index vector at all —
  # consecutive labels span consecutive coordinates, `ranges` being built by accumulation.
  k1 = get(layout.index, first(labels), 0)
  k1 == 0 && error("JointLayout has no coordinates for $(first(labels))")
  consecutive = true
  for (i, label) in enumerate(labels)
    k = get(layout.index, label, 0)
    k == 0 && error("JointLayout has no coordinates for $label")
    if k != k1 + i - 1
      consecutive = false
      break
    end
  end
  consecutive || return getCoordindices(layout, labels)
  return first(layout.ranges[k1]):last(layout.ranges[k1 + length(labels) - 1])
end

getCoordspan(bundle::DensityBundlePoint, labels::AbstractVector{Symbol}) =
  getCoordspan(bundle.layout, labels)

"""
    $SIGNATURES

Lazily walk the variables, yielding `(; label, statetype, range)` — and `point` as well when there is a
base to read it from.  `labels`, `statetypes` and `ranges` are parallel arrays, and this is the one
place that fact is written down.

A bundle point yields the wider tuple, so a loop written against a layout also runs over a bundle but
not the reverse.  Destructure by name (`v.range`), never by position.
"""
eachvariable(layout::JointLayout) = (
  (; label = layout.labels[k], statetype = getStatetype(layout, k), range = layout.ranges[k]) for
  k in eachindex(layout.labels)
)

eachvariable(bundle::DensityBundlePoint) = (
  (;
    label = bundle.layout.labels[k],
    statetype = getStatetype(bundle.layout, k),
    range = bundle.layout.ranges[k],
    point = bundle.point[k],
  ) for k in eachindex(bundle.layout.labels)
)

"""
    $SIGNATURES
The sub-layout over `labels`, which must be a subset of `layout.labels` **in `layout`'s own order** — see
[`subsetBundle`](@ref).
"""
subsetLayout(layout::JointLayout, labels::AbstractVector{Symbol}) = JointLayout(
  collect(Symbol, labels),
  DataType[getStatetype(layout, layout.index[label]) for label in labels],
)

"""
    $SIGNATURES

`labels` re-blocked into the `ArrayPartition` shape the parametric solvers pass around as
`varlabelsAP`.  The bundle stores them flat because that is the coordinate order.
"""
getLabelpartition(bundle::DensityBundlePoint) = getLabelpartition(bundle.layout)

function getLabelpartition(layout::JointLayout)
  # `blockstarts` is why this needs no running offset — and a running offset here would be a captured,
  # mutated closure variable, which Julia boxes.
  blocks = map(eachindex(layout.blocks)) do b
    start = layout.blockstarts[b]
    layout.labels[start:(start + layout.blocks[b].count - 1)]
  end
  return ArrayPartition(blocks...)
end

"""
    $TYPEDEF

The **conditional** half of an elimination, `p(F│S)`, already solved: `ΔF = v - W·ΔS`, with precision
`Λ` over the frontals that does not depend on `ΔS` at all.

$(TYPEDFIELDS)
"""
struct GaussianConditional
  """frontal precision `Λ_FF`; independent of `ΔS`, and singular on a partial clique"""
  Λ::Matrix{Float64}
  """`Λ_FF⁻¹η_F` — the frontal step taken when `ΔS = 0`."""
  v::Vector{Float64}
  """`Λ_FF⁻¹Λ_FS` — the gain, how the frontals follow the separators"""
  W::Matrix{Float64}
  """which variables this is a density over, in `v`'s coordinate order"""
  frontals::Vector{Symbol}
  """which variables it is conditioned on, in `W`'s column order"""
  separators::Vector{Symbol}
end
