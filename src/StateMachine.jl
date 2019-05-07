
"""
    $TYPEDEF

Generic state machine functor type.

Example
```julia
bar!(usrdata) = IncrementalInference.exitStateMachine
foo!(usrdata) = bar!

sm = StateMachine(next=foo!)
usrdata = nothing
while st(usrdata); end
```

Notes
- Also see IncrementalInference/test/testStateMachine.jl
"""
mutable struct StateMachine{T}
  next::Function
  iter::Int
  history::Vector{Tuple{Int, Function, T}}
  StateMachine{T}(;next=emptyState, iter::Int=0) where T = new{T}(next, iter, Vector{Tuple{Int, Function, T}}())
end

"""
    $SIGNATURES

Run state machine function (as functor).

Example
```julia
bar!(usrdata) = IncrementalInference.exitStateMachine
foo!(usrdata) = bar!

sm = StateMachine(next=foo!)
usrdata = nothing
while st(usrdata); end
```
"""
function (st::StateMachine{T})(userdata::T=nothing;
                               breakafter::Function=exitStateMachine,
                               verbose::Bool=false,
                               iterlimit::Int=-1,
                               recordhistory::Bool=false  ) where {T}
  #
  st.iter += 1
  !verbose ? nothing : println("State machine iter=$(st.iter)")
  retval = st.next != breakafter && (iterlimit == -1 || st.iter < iterlimit)
  recordhistory ? push!(st.history, (st.iter, deepcopy(st.next), deepcopy(userdata))) : nothing
  st.next = st.next(userdata)
  return retval
end

"""
    $SIGNATURES

Dummy function in case a `statemachine.next` is not initialized properly.
"""
function emptyState(dummy)
  @warn "Empty state machine, assign `next` to entry function -- i.e. StateMachine(next=foo)"
  return exitStateMachine
end

"""
    $SIGNATURES

Default function used for exiting any state machine.
"""
function exitStateMachine(dummy)
  return emptyState
end




#
