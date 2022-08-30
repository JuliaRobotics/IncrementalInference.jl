
## TBD, will be redone together with fixes for #1010

"""
    $SIGNATURES

Return a default factor type between two variables of types T1 and T2.

Notes
- Most likely used with deconvolution between target variables
"""
function selectFactorType(
  Modl::Module,
  T1::Type{<:InferenceVariable},
  T2::Type{<:InferenceVariable},
)
  return getfield(Modl, Symbol(T1, T2))
end
function selectFactorType(T1::Type{<:InferenceVariable}, T2::Type{<:InferenceVariable})
  return selectFactorType(typeof(T1()).name.module, T1, T2)
end
selectFactorType(T1::Type{<:Position1}, T2::Type{<:Position1}) = LinearRelative{1}
function selectFactorType(T1::Type{<:Position{N}}, T2::Type{<:Position{N}}) where {N}
  return LinearRelative{N}
end
function selectFactorType(T1::InferenceVariable, T2::InferenceVariable)
  return selectFactorType(typeof(T1), typeof(T2))
end
function selectFactorType(dfg::AbstractDFG, s1::Symbol, s2::Symbol)
  return selectFactorType(getVariableType(dfg, s1), getVariableType(dfg, s2))
end

#
