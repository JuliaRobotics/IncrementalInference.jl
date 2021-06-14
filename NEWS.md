# Additional NEWS on IncrementalInference.jl Releases

Currently general maintenance and bug fix changes are mostly tracked via Github Integrations.  E.g. see Milestones along with Label filters to quickly find specific issues.
- https://github.com/JuliaRobotics/IncrementalInference.jl/milestones?state=closed

Also see automated TagBot Release note, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/releases

Alternatively, either use the Github Blame, or the Github `/compare/v0.18.0...v0.19.0` API, e.g.:
- https://github.com/JuliaRobotics/IncrementalInference.jl/compare/v0.18.0...v0.19.0

The list below highlights major breaking changes, and please note that significant efforts are made to properly deprecate old code/APIs according to normal semver workflow -- i.e. breaking changes go through at least one deprecatation (via warnings) on the dominant number in the version number.  E.g. v0.18 -> v0.19 (warnings) -> v0.20 (breaking).

# Major changes in IIF v0.24

- Update compat for ManifoldsBase.jl v0.11 with `AbstractManifold`.
- Transition to only `getManifold` (instead of `getManifolds`), thereby moving towards exclusively using Manifolds.jl, see #1234.
- Deprecate use of `getFactorMean`, use `IIF.getParametricMeasurement` instead.
- Upstreamed `is/set Marginalized` to DFG (#1269).
# Major changes in IIF v0.23

- New `@defVariable` only uses `ManifoldsBase.Manifold` as base abstraction for variable types.
# Major changes in IIF v0.22

- Work in progress toward `ManifoldsBase.Manifold` as base abstraction for variable types.
# Major changes in IIF v0.21

- `CalcResidual` no longer takes a `residual` as input parameter and should return `residual`, see #467 .

# Major changes in IIF v0.20

- The user factor API call strategy has been simplified via `CalcResidual`, see #467 for details.
- User factor API for `getSample` and `.specialsampler` has been standardized via `CalcResidual` (#927) -- for ongoing work please follow #1099 and #1094 and #1069.