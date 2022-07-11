# Additional NEWS on IncrementalInference.jl Releases

Currently general maintenance and bug fix changes are mostly tracked via Github Integrations.  E.g. see Milestones along with Label filters to quickly find specific issues.
- https://github.com/JuliaRobotics/IncrementalInference.jl/milestones?state=closed

Also see automated TagBot Release note, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/releases

Alternatively, either use the Github Blame, or the Github `/compare/v0.18.0...v0.19.0` API, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/compare/v0.18.0...v0.19.0

The list below highlights breaking changes according to normal semver workflow -- i.e. breaking changes go through at least one deprecatation (via warnings) on the dominant number in the version number.  E.g. v0.18 -> v0.19 (warnings) -> v0.20 (breaking).  Note that ongoing efforts are made to properly deprecate old code/APIs

# Changes in v0.30

- `ArrayPartition` should be used instead of `ProductRepr`, see issue #1537.
- Remove old deprecated option keywords in `addVariable` and `addFactor`.

# Changes in v0.29

- Upgrade to Manifolds.jl v0.8
- Deprecate `initManual!`, instead use `initVariable!`.
# Changes in v0.28

- `HeatmapGridDensity` now only supports `ManifoldKernelDensity` functions.
- `PackedHeatmapGridDensity` has an expanded fields to support future stash and cache serialization strategies.
- Internal `parchDistribution` functions have been added towards future stashed serialization strategies.
- Internal `_update!` function supports updating of the `HeatmapGridDensity` distribution.
- Unpacking of `PackedManifoldKernelDensity` is more versatile with improved `.partial` and `.bw` options.
- Bugfix on `multihypo=` which now includes `nullSurplus` on sibling relative factors to a variable with a `multihypo` factor, #1518.

# Changes in v0.27

- InMemDFGType is deprecated in favor of LocalDFG (exported from DistributedFactorGraphs).
- Factor serialization is now top level JSON only #1476.
- Serialization of distributions are now JSON only #1468, #1472, #1473 (removed custom string legacy).
- Fix chicken and egg problem on unpackFactor, change `convert` to `reconstFactorData`, #1424.
- Add factor `preambleCache(dfg, vecVars, usrfnc)`, #1462, #1466.  Doesn't work for parametric yet (#1480).
- Add `CalcFactor.cache` using preamble, #1481.  Not thread safe yet.
- Standardize local graph naming to `LocalDFG`, #1479.
- Refactor getDimension and sampling, #1463.
- Language upgrades on `qr` for Julia 1.7, #1464.
- Various other fixes and upgrades, https://github.com/JuliaRobotics/IncrementalInference.jl/milestone/111?closed=1
- Add distribution serialization for Rayleigh.
- Add `Position{N}` and `Position1`..`Position4` as new standard and aliases for `ContinuousScalar`, `ContinuousEuclid{N}`.

# Changes in v0.26

- Standarding (non-binding) easy factor dipatch cases so measurement field is under `.Z` (#1441).
- `CalcFactor._allowThreads` can now be used as workaround for `Threads` yield blocking issue during first run (#1451).
- Canonical graph generator API change to `generateGraph_ABC` (#1454).

# Changes in v0.25

- Changed API to `testFactorResidualBinary(fct, meas::Tuple, (T_i, param_i),...)` to grow beyond binary.
- PPE methods used keyword `method::AbstractPointParametricType` which is now replaced with the keyword `ppeType`.
- Belief points are now stored as a `Vector{P}` (rather than legacy `Matrix`), and currently still under the restriction `P <: AbstractVector{<:Real}`.  Objective is moving to `P` any user desired point that fits with the JuliaManifolds/Manifolds.jl patterns.
- Deprecating use of `ensureAllInitialized!`, use `initAll!` instead.
- Upstream `calcHelix_T` canonical generator utility from RoME.jl.
- Deserialization of factors with DFG needs new API and change of solverData and CCW type in factor.
- Deprecate use of `getParametricMeasurement` and use `getMeasurementParametric` instead, and add `<:AbstractManifold` to API.
- Deprecate use of `solveBinaryFactorParameteric`, instead use `solveFactorParameteric`.
- Deprecating `approxConvBinary`, use `approxConvBelief` instead.
- Removing obsolete `approxConvCircular`, use `approxConvBelief` instead.
- `getSample` should return a single sample and no longer takes the N(number of samples) parameter.
- `solveTree!` / `solveGraph!` now returns just one value `tree<:AbstractBayesTree`.  Previous version returned three values, `tree, smt, hist` (#1379).
- **Note for v0.25.5** Serialization of newly introduced type `PackedHeatmapGridDensity` changed from v0.25.4, unlikely have yet been used publically, therefore emphasizing fastest possible standardization in this case (even though this particular event does not strictly follow semver).  General usage and operation is effectively unchanged,see #1435.

# Changes in v0.24

- Update compat for ManifoldsBase.jl v0.11 with `AbstractManifold`.
- Transition to only `getManifold` (instead of `getManifolds`), thereby moving towards exclusively using Manifolds.jl, see #1234.
- Deprecate use of `getFactorMean`, use `IIF.getParametricMeasurement` instead.
- Upstreamed `is/set Marginalized` to DFG (#1269).
# Changes in v0.23

- New `@defVariable` only uses `ManifoldsBase.Manifold` as base abstraction for variable types.
# Changes in v0.22

- Work in progress toward `ManifoldsBase.Manifold` as base abstraction for variable types.
# Changes in v0.21

- `CalcResidual` no longer takes a `residual` as input parameter and should return `residual`, see #467 .

# Changes in v0.20

- The user factor API call strategy has been simplified via `CalcResidual`, see #467 for details.
- User factor API for `getSample` and `.specialsampler` has been standardized via `CalcResidual` (#927) -- for ongoing work please follow #1099 and #1094 and #1069.
