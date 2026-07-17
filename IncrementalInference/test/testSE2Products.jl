#=
Reference & regression test: fusing (multiplying) several Gaussian priors on a 2D pose.

This file compares three independent ways of computing the product of N Gaussian
priors placed on a 2D pose manifold, and cross-checks them against each other:

  - `prod_by_fg`         : build a 1-variable factor graph with N `ManifoldPrior` factors
                            and let the parametric solver fuse them.
  - `prod_by_bruteforce`  : evaluate the product density on a dense grid, take its mode,
                            then fit a Gaussian to the (weighted) grid samples. Slow, but
                            makes no linearization assumptions -- treated as ground truth.
  - `prod_by_amp`         : the closed-form on-manifold Gaussian product used internally
                            by ApproxManifoldProducts (`calcProductGaussians`).

It also compares two different *manifolds* for representing a 2D pose:

  - `Pose2_SE2` : the true SE(2) Lie group (semidirect product, right-variant).
  - `Pose2`     : the product manifold `TranslationGroup(2) × SpecialOrthogonalGroup(2)`

`check_mean_invariance` makes the difference concrete: it perturbs a set of poses by a
fixed rigid transform `h` (left- and right-multiplication, and inversion) and checks
whether "fuse-then-transform" equals "transform-then-fuse". For `Pose2_SE2`, SE(2) has
no bi-invariant metric, so this only holds exactly on the side matching the group's
trivialization. For `Pose2`, `compose` is not a real rigid-body transform at all, so the check
holds trivially on both sides -- that "pass" is not a sign of correctness, see the
testset below for details.
=#

using Test
using IncrementalInference
using DistributedFactorGraphs
using Distributions
using LieGroups
using LinearAlgebra
using StaticArrays
import ApproxManifoldProducts
import Rotations as _Rot

##==============================================================================
## Two state types representing a 2D pose, for comparison
##==============================================================================

# True SE(2) Lie group: translation and rotation are coupled through `compose`.
DFG.@defStateType(
    Pose2_SE2,
    SpecialEuclideanGroup(2; variant = :right),
    ArrayPartition(@SVector([0.0, 0.0]), @SMatrix([1.0 0.0; 0.0 1.0])),
)

# Product manifold: translation and rotation are independent under `compose`.
# (This mirrors how `RoME.Pose2` is defined -- kept local here so this file has no
# dependency on RoME.)
DFG.@defStateType(
    Pose2,
    TranslationGroup(2) × SpecialOrthogonalGroup(2),
    ArrayPartition(@SVector([0.0, 0.0]), @SMatrix([1.0 0.0; 0.0 1.0])),
)

##==============================================================================
## Product-of-Gaussians via each method
##==============================================================================

"""
    prod_by_fg(statekind, points, covars)

Fuse Gaussian priors `(points[i], covars[i])` by building a 1-variable factor graph
(one `ManifoldPrior` per input) and running the parametric solver. Returns `(μ, Σ)`.
"""
function prod_by_fg(statekind, points, covars)
    length(points) == length(covars) || error("points and covars must have the same length")
    !isempty(points) || error("points and covars must be non-empty")

    fg = initfg()
    fg.solverParams.graphinit = false
    addVariable!(fg, :x0, statekind)

    G = getManifold(fg, :x0)
    for (pt, Σ) in zip(points, covars)
        addFactor!(fg, [:x0], ManifoldPrior(G, pt, MvNormal(Σ)))
    end

    IIF.solveGraphParametric!(
        fg;
        init = true,
        is_sparse = false,
        finiteDiffCovariance = true,
        jacobian_method = :forwarddiff,
        damping_term_min = 1e-3,
    )

    μ = getBelief(fg, :x0, :parametric).points[1]
    Σ_μ = getBelief(fg, :x0, :parametric).trailing_forms[1]

    return μ, Σ_μ
end

"""
    prod_by_amp(statekind, points, covars)

Fuse Gaussian priors using ApproxManifoldProducts' closed-form on-manifold Gaussian
product (`calcProductGaussians`), i.e. no iterative solve. Returns `(μ, Σ)`.
"""
function prod_by_amp(statekind, points, covars)
    G = getManifold(statekind)
    m, Σ = ApproxManifoldProducts.calcProductGaussians(G, points, covars)
    return m, Σ
end

"""
    prod_by_bruteforce(statekind, points, covars; xs, ys, θs, x_step, y_step, θ_step)

Ground-truth-ish product: evaluate the product of the input densities on a dense
(x, y, θ) grid, take the mode, then fit a Gaussian (in the tangent space at the mode)
to the weighted grid samples, re-linearizing once at the fitted mean. No search-window
kwarg needs to be given -- defaults are estimated from the input points/covariances --
but pass `xs`/`ys`/`θs` explicitly to keep the grid (and therefore runtime) small.
Returns `(μ, Σ, details)` where `details` exposes the grid and per-point densities.
"""
function prod_by_bruteforce(statekind, points, covars; xs = nothing, ys = nothing, θs = nothing, x_step = 0.1, y_step = 0.1, θ_step = 0.01)

    length(points) == length(covars) || error("points and covars must have the same length")
    !isempty(points) || error("points and covars must be non-empty")

    G = getManifold(statekind)
    lieG = LieAlgebra(G)

    # Estimate default x/y search windows from point locations and covariance spread.
    tx = map(pt -> pt.x[1][1], points)
    ty = map(pt -> pt.x[1][2], points)
    x_center = sum(tx) / length(tx)
    y_center = sum(ty) / length(ty)
    σx_max = maximum(map(Σ -> sqrt(max(Σ[1, 1], eps(Float64))), covars))
    σy_max = maximum(map(Σ -> sqrt(max(Σ[2, 2], eps(Float64))), covars))
    x_span = maximum(tx) - minimum(tx)
    y_span = maximum(ty) - minimum(ty)

    if isnothing(xs)
        halfspan_x = max(3.0, 0.5 * x_span + 4.0 * σx_max)
        xs = (x_center - halfspan_x):x_step:(x_center + halfspan_x)
    end
    if isnothing(ys)
        halfspan_y = max(3.0, 0.5 * y_span + 4.0 * σy_max)
        ys = (y_center - halfspan_y):y_step:(y_center + halfspan_y)
    end

    # Center the default theta search window around the circular mean of input headings.
    if isnothing(θs)
        angles = map(points) do pt
            R = pt.x[2]
            atan(R[2, 1], R[1, 1])
        end
        θ_center = atan(sum(sin, angles), sum(cos, angles))

        # Expand theta support using both heading spread and covariance scale.
        angle_offsets = map(a -> atan(sin(a - θ_center), cos(a - θ_center)), angles)
        θ_spread = maximum(abs, angle_offsets)
        σθ_max = maximum(map(Σ -> sqrt(max(Σ[3, 3], eps(Float64))), covars))
        θ_halfspan = min(pi, max(1.0, θ_spread + 5.0 * σθ_max))

        θs = (θ_center - θ_halfspan):θ_step:(θ_center + θ_halfspan)
    end

    grid_points = map(Iterators.product(xs, ys, θs)) do (x, y, θ)
        ArrayPartition(SA[x, y], _Rot.RotMatrix2(θ))
    end

    pdf_terms = map(zip(points, covars)) do (pt, Σ)
        map(grid_points) do gp
            X = log(G, pt, gp)
            Xc_e = vee(lieG, X)
            pdf(MvNormal(Σ), Xc_e)
        end
    end

    pdf_prod = ones(size(grid_points))
    for pdf_i in pdf_terms
        pdf_prod .*= pdf_i
    end

    _, mode_index = findmax(pdf_prod)
    mode_point = grid_points[mode_index]

    tangent_samples = map(grid_points) do gp
        vee(lieG, log(G, mode_point, gp))
    end
    X = reduce(hcat, tangent_samples)
    w = (pdf_prod ./ sum(pdf_prod))[:]
    fit_mvn = fit_mle(MvNormal, X, w)

    δμ = fit_mvn.μ
    μ0 = exp(G, mode_point, hat(lieG, SA[δμ...], ArrayPartition))

    # Re-linearize at μ0 to get a more reliable covariance, especially on theta.
    tangent_samples_μ0 = map(grid_points) do gp
        vee(lieG, log(G, μ0, gp))
    end
    X_μ0 = reduce(hcat, tangent_samples_μ0)
    fit_mvn_μ0 = fit_mle(MvNormal, X_μ0, w)

    δμ_μ0 = fit_mvn_μ0.μ
    μ = exp(G, μ0, hat(lieG, SA[δμ_μ0...], ArrayPartition))
    Σ_μ = Matrix(cov(fit_mvn_μ0))

    details = (; xs, ys, θs, grid_points, pdf_terms, pdf_prod, mode_index, mode_point, fit_mvn, fit_mvn_μ0)
    return μ, Σ_μ, details
end

##==============================================================================
## Equivariance check: does fusion commute with a rigid-body transform?
##==============================================================================

"""
    check_mean_invariance(G, points, covars, h, prod_method; atol, label)

`prod_method(points, covars) -> μ` (or `(μ, Σ, ...)`, only `μ` is used) fuses a set of
poses to a single mean. Checks that fusing commutes with composing every input by a
fixed transform `h`, both on the left (`h * pt`) and right (`pt * h`), and with
inversion (`inv(pt)`). Prints a pass/fail per check and returns a `NamedTuple`.
"""
function check_mean_invariance(G, points, covars, h, prod_method; atol = 1e-6, label = "")

    _mean(points_, covars_) = begin
        out = prod_method(points_, covars_)
        out isa Tuple ? out[1] : out
    end

    m = _mean(points, covars)

    points_left = map(pt -> compose(G, h, pt), points)
    m_left_fused = _mean(points_left, covars)
    m_left_expected = compose(G, h, m)
    left_check = isapprox(G, m_left_fused, m_left_expected; atol)

    points_right = map(pt -> compose(G, pt, h), points)
    m_right_fused = _mean(points_right, covars)
    m_right_expected = compose(G, m, h)
    right_check = isapprox(G, m_right_fused, m_right_expected; atol)

    points_inv = map(pt -> inv(G, pt), points)
    m_inv_fused = _mean(points_inv, covars)
    m_inv_expected = inv(G, m)
    inv_check = isapprox(G, m_inv_fused, m_inv_expected; atol)

    if !isempty(label)
        println("--- ", label, " ---")
    end
    println("Left-Invariance Pass? -> ", left_check)
    println("Right-Invariance Pass? -> ", right_check)
    println("Inverse-Invariance Pass? -> ", inv_check)

    return (; mean = m, left = left_check, right = right_check, inverse = inv_check)
end

##==============================================================================
## Tests
##==============================================================================

@testset "SE(2) pose fusion: FG solve vs brute-force vs AMP closed-form" begin
    p = ArrayPartition(SA[1.0, 2.0], _Rot.RotMatrix2(0.3))
    Σp = diagm(SA[0.20, 0.20, 0.05] .^ 2)

    q = ArrayPartition(SA[1.4, 1.6], _Rot.RotMatrix2(0.6))
    R2 = _Rot.RotMatrix2(pi / 6)
    Rblock = [R2 zeros(2, 1); zeros(1, 2) 1.0]
    Σq = Rblock * diagm([0.15, 0.30, 0.07] .^ 2) * Rblock'

    xs = -1.0:0.05:3.5
    ys = 0.5:0.05:3.5
    θs = 0.0:0.01:1.0

    @testset "Pose2_SE2 (coupled SE(2) Lie group)" begin
        G = getManifold(Pose2_SE2)
        μ_fg, Σ_fg = prod_by_fg(Pose2_SE2, [p, q], [Σp, Σq])
        μ_amp, Σ_amp = prod_by_amp(Pose2_SE2, [p, q], [Σp, Σq])
        μ_bf, Σ_bf, _ = prod_by_bruteforce(Pose2_SE2, [p, q], [Σp, Σq]; xs, ys, θs)

        @test isapprox(G, μ_fg, μ_bf; atol = 5e-2)
        @test isapprox(G, μ_amp, μ_bf; atol = 5e-2)
        @test isapprox(Σ_fg, Σ_bf; atol = 5e-2)
        @test isapprox(Σ_amp, Σ_bf; atol = 5e-2)
    end

    @testset "Pose2 (decoupled product manifold)" begin
        G = getManifold(Pose2)
        μ_fg, Σ_fg = prod_by_fg(Pose2, [p, q], [Σp, Σq])
        #TODO needs AMP v0.15.7
        # μ_amp, Σ_amp = prod_by_amp(Pose2, [p, q], [Σp, Σq])
        μ_bf, Σ_bf, _ = prod_by_bruteforce(Pose2, [p, q], [Σp, Σq]; xs, ys, θs)

        @test isapprox(G, μ_fg, μ_bf; atol = 5e-2)
        @test isapprox(Σ_fg, Σ_bf; atol = 5e-2)
        @test_broken isapprox(G, μ_amp, μ_bf; atol = 5e-2)
        @test_broken isapprox(Σ_amp, Σ_bf; atol = 5e-2)
    end
end

@testset "SE(2) pose fusion: equivariance under rigid-body transform" begin
    p = ArrayPartition(SA[10.0, 18.0], _Rot.RotMatrix2(0.0))
    Σp = diagm(SA[1.0, 1.0, 1.0] .^ 2)

    q = ArrayPartition(SA[10.0, 24.0], _Rot.RotMatrix2(1.0))
    Σq = diagm(SA[1.0, 1.0, 1.0] .^ 2)

    h = ArrayPartition(SA[2.5, -4.0], _Rot.RotMatrix2(0.5))

    @testset "Pose2_SE2 (coupled SE(2) group, :right variant): left-equivariant only" begin
        G = getManifold(Pose2_SE2)
        fg_mean_method = (pts, Σs) -> prod_by_fg(Pose2_SE2, pts, Σs)
        r = check_mean_invariance(G, [q, p], [Σq, Σp], h, fg_mean_method; label = "Pose2_SE2 FG mean")
        # SE(2) has no bi-invariant metric, so a Gaussian-on-manifold fusion can only be
        # equivariant on the side that matches the group's chosen trivialization. `Pose2_SE2`
        # LieGroups.jl is left-trivialized.
        @test r.left
        @test !r.right
        @test !r.inverse
    end

    @testset "Pose2 (decoupled product manifold): equivariant on both sides, trivially" begin
        G = getManifold(Pose2)
        fg_mean_method = (pts, Σs) -> prod_by_fg(Pose2, pts, Σs)
        r = check_mean_invariance(G, [q, p], [Σq, Σp], h, fg_mean_method; label = "Pose2 FG mean")
        # `compose` on the product manifold adds translations and composes rotations independently
        @test r.left
        @test r.right
        @test r.inverse
    end
end

##==============================================================================
## Reference (not run): visualize the brute-force product density and geodesics
##==============================================================================
# Paste this block into a Makie-loaded REPL (`using GLMakie` or `using CairoMakie`)
# to see what `prod_by_bruteforce` is integrating over, and how the Lie-group geodesic
# between `p` and `q` compares to the base-manifold (Riemannian) geodesic.
if false
    using GLMakie

    p = ArrayPartition(SA[10.0, 18.0], _Rot.RotMatrix2(0.0))
    Σp = diagm(SA[3.0, 3.0, 1.0] .^ 2)
    q = ArrayPartition(SA[10.0, 24.0], _Rot.RotMatrix2(1.0))
    Σq = diagm(SA[3.0, 3.0, 1.0] .^ 2)

    G = getManifold(Pose2_SE2)
    _, _, bf_details = prod_by_bruteforce(Pose2_SE2, [p, q], [Σp, Σq])
    xs, ys, θs = bf_details.xs, bf_details.ys, bf_details.θs
    pdf_ps, pdf_qs = bf_details.pdf_terms
    pdf_pqs = bf_details.pdf_prod # brute-force product

    cols = Makie.wong_colors()
    contour(xs, ys, sum(pdf_pqs; dims=3)[:, :, 1]; color = cols[1], linewidth = 2, levels = 15, axis =(aspect=DataAspect(),))
    contour!(xs, ys, sum(pdf_ps; dims=3)[:, :, 1]; color = cols[2], alpha = 0.5)
    contour!(xs, ys, sum(pdf_qs; dims=3)[:, :, 1]; color = cols[3], alpha = 0.5)

    # Lie-group geodesic between p and q (on G)
    t_geo = range(0.0, 1.0; length = 100)
    pq_geo = map(t_geo) do t
        exp(G, p, t * log(G, p, q))
    end
    pq_geo_x = [first(g.x)[1] for g in pq_geo]
    pq_geo_y = [first(g.x)[2] for g in pq_geo]
    lines!(pq_geo_x, pq_geo_y; color = cols[4], linewidth = 3) # magenta

    # untested to plot the Riemannian geodesic (straight line), we use the base manifold
    M = base_manifold(G)
    t_geo = range(0.0, 1.0; length = 100)
    pq_geo = map(t_geo) do t
        exp(M, p, t * log(M, p, q))
    end
    pq_geo_x = [first(g.x)[1] for g in pq_geo]
    pq_geo_y = [first(g.x)[2] for g in pq_geo]
    lines!(pq_geo_x, pq_geo_y; color = cols[5], linewidth = 3)
end
