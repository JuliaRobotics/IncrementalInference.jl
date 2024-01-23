# NEWS on IncrementalInference.jl Releases

Currently general maintenance and bug fix changes are mostly tracked via Github Integrations.  E.g. see Milestones along with Label filters to quickly find specific issues.
- https://github.com/JuliaRobotics/IncrementalInference.jl/milestones?state=closed

Also see automated TagBot Release note, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/releases

Alternatively, either use the Github Blame, or the Github `/compare/v0.18.0...v0.19.0` API, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/compare/v0.18.0...v0.19.0

The list below highlights breaking changes according to normal semver workflow -- i.e. breaking changes go through at least one deprecatation (via warnings) on the dominant number in the version number.  E.g. v0.18 -> v0.19 (warnings) -> v0.20 (breaking).  Note that ongoing efforts are made to properly deprecate old code/APIs

# Changes in v0.35

- Standardize toward Manopt.jl (currently Riemannian Levenberg-Marquart), still have Optim.jl legacy support (#1784, #1778).
- Much faster solves, both parametric and nonparametric (#1783, #1782, #1793).
- Better standardize relative factors to use of tangent vectors (#1790).
- Now abstract type `CalcFactor` with dedicated dispatches in various cases, e.g. `CalcFactorNormSq`, etc. (#1786).
- Bug fixes and inference improvements (#1781, #1785, #1789)
- Support for Julia 1.10.
- Extension usage of AMD.jl for `ccolamd` variable ordering features, dropped internal SuiteSparse calls.  JL 1.10 removal of `SuiteSparse_long` (#1763).
- Further bug fixes for transition to `StaticArrays` value stores and computes, including `Position{N}` (#1779, #1776).
- Restore `DifferentialEquation.jl` factor `DERelative` functionality and tests that were suppressed in a previous upgrade (#1774, #1777).
- Restore previously suppressed tests (#1781, #1721, #1780)
- Improve DERelative factor on-manifold operations (#1775, #1802, #1803).
- Fixed a typo via deprecation, `solveFactorParametric` replaces `solveFactorParameteric`.

# Changes in v0.34

- Start transition to Manopt.jl via Riemannian Levenberg-Marquart.
- Deprecate `AbstractRelativeRoots`.
- Standardization improvements surrounding weakdeps code extensions. 
- Code quality improvements along wiht refactoring and reorganizing of file names and locations.
- Restoring `DERelative` factors, through v0.34.1 and v0.34.2.
- Switching to weakdep AMD.jl for `ccolmod` dependency, part of Julia 1.10 upgrade.  Dropping `SuiteSparse_long` dependency.  Further fixes necessary to restore full user constrained tree variable order functionality.

# Changes in v0.33

- Upgrades for DFG using StructTypes.jl (serialization).
# Changes in v0.32

- Major internal refactoring of `CommonConvWrapper` to avoid abstract field types, and better standardization; towards cleanup of internal multihypo handling and naming conventions.
- Internal refactoring removing several legacy fields from `CalcFactor`.
- All factors now require a `getManifold` definition.
- Now have `CalcFactor.manifold` to reduce new allocation load inside hot-loop for solving.
- Fixed tesing issues in `testSpecialEuclidean2Mani.jl`.
- Refactored, consolidated, and added more in-place operations in surrounding `ccw.measurement`.
- Refactor `CommonConvWrapper` to a not non-mutable struct, with several cleanups and some updated compat requirements.
- Refactor interal hard type `HypoRecipe`.
- Add `MetaPrior` for adding meta data but not influencing the numerical solution.

# Changes in v0.31
- `FactorMetaData` is deprecated and replaced by `CalcFactor`.
- Updated `Base.deepcopy_internal` fix for use with Julia 1.8.1, see #1629.
- Added a few missing exports incl. `getTags`, `_update!, see #1626 #1628.
- Refactoring to remove `FactorMetadata` (#1611) and `ConvPerThread` (#1615, #1625) objects, which is consolidated into `CalcFactor` and `CommonConvWrapper`.
- Added JuliaFormatter, see #1620.
- Add `SnoopPrecompile.jl` on a few basic solve features to start, see #1631.
- Support n-ary parametric solving such as OAS factors.

# Changes in v0.30

- `ArrayPartition` should be used instead of `ProductRepr`, see issue #1537.
- Remove old deprecated option keywords in `addVariable` and `addFactor`.
- Improve `IIF.solveGraphParametric`.
- Introduce `IIF.autoinitParametric!`.
- Upgrade `initAll!(dfg, :parametric)`.
- Refactor many files to subfolders `src/services` or `src/entities`.

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
- Deprecate use of `solveBinaryFactorParameteric`, instead use `solveFactorParametric`.
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
