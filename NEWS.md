# Additional NEWS on IncrementalInference.jl Releases

Currently general maintenance and bug fix changes are mostly tracked via Github Integrations.  E.g. see Milestones along with Label filters to quickly find specific issues.
- https://github.com/JuliaRobotics/IncrementalInference.jl/milestones?state=closed

Also see automated TagBot Release note, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/releases

Alternatively, either use the Github Blame, or the Github `/compare/v0.18.0...v0.19.0` API, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/compare/v0.18.0...v0.19.0

The list below highlights major breaking changes, and please note that significant efforts are made to properly deprecate old code/APIs according to normal semver workflow -- i.e. breaking changes go through at least one deprecatation (via warnings) on the dominant number in the version number.  E.g. v0.18 -> v0.19 (warnings) -> v0.20 (breaking).

# Major changes in v0.25

- Changed API to `testFactorResidualBinary(fct, meas::Tuple, (T_i, param_i),...)` to grow beyond binary.
- PPE methods used keyword `method::AbstractPointParametricType` which is now replaced with the keyword `ppeType`.
- Belief points are now stored as a `Vector{P}` (rather than legacy `Matrix`), and currently still under the restriction `P <: AbstractVector{<:Real}`.  Objective is moving to `P` any user desired point that fits with the JuliaManifolds/Manifolds.jl patterns.
- Deprecating use of `ensureAllInitialized!`, use `initAll!` instead.
- New helper function `randToPoints(::SamplableBelief, N=1)::Vector{P}` to help with `getSample` for cases with new `ManifoldKernelDensity` beliefs for manifolds containing points of type `P`.
- Upstream `calcHelix_T` canonical generator utility from RoME.jl.
- Deserialization of factors with DFG needs new API and change of solverData and CCW type in factor.
- Deprecate use of `getParametricMeasurement` and use `getMeasurementParametric` instead, and add `<:AbstractManifold` to API.
- Deprecate use of `solveBinaryFactorParameteric`, instead use `solveFactorParameteric`.
- Deprecating `approxConvBinary`, use `approxConvBelief` instead.
- Deprecating `accumulateFactorChain`, use `approxConvBelief` instead.


# Major changes in v0.24

- Update compat for ManifoldsBase.jl v0.11 with `AbstractManifold`.
- Transition to only `getManifold` (instead of `getManifolds`), thereby moving towards exclusively using Manifolds.jl, see #1234.
- Deprecate use of `getFactorMean`, use `IIF.getParametricMeasurement` instead.
- Upstreamed `is/set Marginalized` to DFG (#1269).
# Major changes in v0.23

- New `@defVariable` only uses `ManifoldsBase.Manifold` as base abstraction for variable types.
# Major changes in v0.22

- Work in progress toward `ManifoldsBase.Manifold` as base abstraction for variable types.
# Major changes in v0.21

- `CalcResidual` no longer takes a `residual` as input parameter and should return `residual`, see #467 .

# Major changes in v0.20

- The user factor API call strategy has been simplified via `CalcResidual`, see #467 for details.
- User factor API for `getSample` and `.specialsampler` has been standardized via `CalcResidual` (#927) -- for ongoing work please follow #1099 and #1094 and #1069.