
## TBD, will be redone together with fixes for #1010

"""
    $SIGNATURES

Return a default factor type between two variables of types T1 and T2.

Notes
- Most likely used with deconvolution between target variables
"""
selectFactorType(Modl::Module, T1::Type{<:InferenceVariable}, T2::Type{<:InferenceVariable}) = getfield(Modl, Symbol(T1, T2))
selectFactorType(T1::Type{<:InferenceVariable}, T2::Type{<:InferenceVariable}) = selectFactorType(typeof(T1()).name.module, T1, T2)
selectFactorType(T1::Type{<:Position1}, T2::Type{<:Position1}) = LinearRelative{1}
selectFactorType(T1::Type{<:Position{N}}, T2::Type{<:Position{N}}) where N = LinearRelative{N}
selectFactorType(T1::InferenceVariable, T2::InferenceVariable) = selectFactorType(typeof(T1), typeof(T2))
selectFactorType(dfg::AbstractDFG, s1::Symbol, s2::Symbol) = selectFactorType( getVariableType(dfg, s1), getVariableType(dfg, s2) )

#