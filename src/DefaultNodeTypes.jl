
## TBD, will be redone together with fixes for #1010

"""
    $SIGNATURES

Hacky version to return which factor type to use between two variables of types T1 and T2.

Pending work on RoME.jl #244 and IIF.jl #1010 and others, after which this will be refactored.
"""
selectFactorType(T1::InstanceType{ContinuousScalar}, T2::InstanceType{ContinuousScalar}) = LinearRelative{1}
selectFactorType(T1::Type{<:ContinuousEuclid{N}}, T2::Type{<:ContinuousEuclid{N}}) where N = LinearRelative{N}
selectFactorType(T1::InferenceVariable, T2::InferenceVariable) = selectFactorType(typeof(T1), typeof(T2))
selectFactorType(dfg::AbstractDFG, s1::Symbol, s2::Symbol) = selectFactorType( getVariableType(dfg, s1), getVariableType(dfg, s2) )

