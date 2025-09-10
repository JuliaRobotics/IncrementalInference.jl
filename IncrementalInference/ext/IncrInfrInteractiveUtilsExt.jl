module IncrInfrInteractiveUtilsExt

@info "IncrementalInference.jl is loading extension related to InteractiveUtils.jl."

using InteractiveUtils
using DocStringExtensions
using IncrementalInference: StateType, AbstractPriorObservation, RelativeObservation
# using IncrementalInference: getCurrentWorkspaceFactors, getCurrentWorkspaceVariables, listTypeTree
import IncrementalInference: getCurrentWorkspaceFactors, getCurrentWorkspaceVariables, listTypeTree

export getCurrentWorkspaceFactors, getCurrentWorkspaceVariables
export listTypeTree

"""
    $(SIGNATURES)

Return all factors currently registered in the workspace.
"""
function getCurrentWorkspaceFactors()
  return [
    InteractiveUtils.subtypes(AbstractPriorObservation)...,
    InteractiveUtils.subtypes(RelativeObservation)...,
  ]
end

"""
    $(SIGNATURES)

Return all variables currently registered in the workspace.
"""
function getCurrentWorkspaceVariables()
  return InteractiveUtils.subtypes(StateType)
end

function _listTypeTree(mytype, printlevel::Int)
  allsubtypes = InteractiveUtils.subtypes(mytype)
  for cursubtype in allsubtypes
    print("\t"^printlevel)
    println("|___", cursubtype)
    printlevel += 1
    _listTypeTree(cursubtype, printlevel)
    printlevel -= 1
  end
end

"""
    $SIGNATURES
List the types that inherit from `T`.

Notes
- from https://youtu.be/S5R8zXJOsUQ?t=1531
"""
function listTypeTree(T)
  println(T)
  return _listTypeTree(T, 0)
end

end #module