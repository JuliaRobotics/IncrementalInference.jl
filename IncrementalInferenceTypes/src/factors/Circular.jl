
"""
$(TYPEDEF)

Factor between two Sphere1 variables.

Related

[`Sphere1`](@ref), [`PriorSphere1`](@ref), [`Polar`](@ref), [`ContinuousEuclid`](@ref)
"""
DFG.@defObservationType CircularCircular RelativeObservation LieGroups.CircleGroup()


"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Circular variable:

Example:
--------
```julia
PriorCircular( MvNormal([10; 10; pi/6.0], diagm([0.1;0.1;0.05].^2)) )
```

Related

[`Circular`](@ref), [`Prior`](@ref), [`PartialPrior`](@ref)
"""
DFG.@defObservationType PriorCircular AbstractPrior LieGroups.CircleGroup()
