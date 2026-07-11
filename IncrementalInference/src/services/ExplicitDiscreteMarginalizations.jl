# Second iteration of explicitly exloring the marginalization of discrete variables onto the continuous state space.
# Although this code is still excessive and messy, this is a significant feature expansion; and future versions will
# generalize the marginalization process to allow for implicit hypothesis exploration.  The messy explicit version of code is
# intended to help develop the required generalistic unit tests.  The unit tests will then be validated, frozen and used to
# confirm that future "algebraic" marginalization (implicit) versions operate correctly.
# FYI, the complexity of general multihypothesis convolutions can be deceiving, however, note that the coding
# complexity is contained for each indivual factor at a time.  Global Bayes tree inference then creates the symphony
# of non-Gaussian (multimodal) posterior beliefs from the entire factor graph.
#
# 2018/6/01 @dehann

"""
    $SIGNATURES

Return common vectors `(allmhp, certainidx,uncertnidx)` used for dealing with multihypo cases.
"""
function getHypothesesVectors(
  mhp::Vector{Float64},
)::Tuple{Vector{Int}, Vector{Int}, Vector{Int}}
  allmhp = 1:length(mhp)
  certainidx = allmhp[mhp .== 0.0]  # TODO remove after gwp removed
  uncertnidx = allmhp[0.0 .< mhp]
  return (allmhp, certainidx, uncertnidx)
end



"""
    $(SIGNATURES)

This function explicitly encodes the marginalization of a discrete categorical selection variable 
for ambiguous data association situations.  This function populates `allelements` with particle 
indices associated with particular multihypothesis selection while `activehypo` simultaneously 
contains the hypothesis index and factor graph variables associated with that hypothesis selection.  The return value `certainidx` are the hypotheses that are not in question.

This function does not consider whether a unroll hypo lambda can actually exist, it just generates
a recipe of which options should be considered.  Whoever solves the recipe is responsible for 
building the right lambda or doing something else.

Input:
- `maxlen` is the max number of samples across all variables

Output:
- `certainidx`:   non fractional variables
- `allelements`:  list of which particles go with which hypothesis selection
- `activehypo`:   list of which hypotheses go with which certain + fractional variables
- `mhidx`:        multihypothesis selection per particle idx
- `sfidx`:        Solve for idx

Example:
```julia
idx=(1,2,3)
multihypo=[1.0;0.5;0.5]    # X,La,Lb

# specfic example
# this is important -- e.g. `pts = approxConv(fg, :XLaLbf1, :La)`
sfidx=2

certainidx=
1-element Array{Int64,1}:
 1                         # X
allelements=
3-element Array{Any,1}:
 Int64[]
 [1, 2, 11, ...]
 [3, 4, 5, ...]
activehypo=
3-element Array{Any,1}:
 (0, [2])     # nullhypo -- forced afterward, might be deprecated if better solution is found
 (1, [1, 2])  # unroll hypo lambda for X,La
 (2, [1, 2])  # unroll hypo lambda for X,La
 (3, [2, 3])  # would be (but cannot) unroll hypo La,Lb # almost our nullhypo # wont build a lambda

# now select which lambdas to build based on converted `activehypo`  rand(  Categorical( [0; 0.5;0.5] ) )
 mhidx=
100-element Array{Int64,1}:
2, 2, 3, 3, 3,...
```

Another example on what is `certainidx`
```julia
multihypo=[1;1;0.5;0.5]
# results in: `certainidx = [1;2]`
```

Notes:
- Issue 427, race condition during initialization since n-ary variables not resolvable without other init.

DevNotes
- FIXME convert into some kind of `HypoRecipe` struct, to improve readibility and code maintainability
- TODO add nullhypo cases to returning result
- FIXME make type-stable `activehypo` and others
- Improved implementations should implicitly induce the same behaviour through summation (integration) when marginalizing any number of discrete variables.

```
# `allelements` example BearingRange [:x1, 0.5:l1a, 0.5:l1b]
# sfidx = (1=:x1,2=:l1a,3=:l1b)
if solvefor :x1, then allelem = [mhidx.==:l1a; mhidx.==l1b]
if solvefor :l1a, then allelem = [mhidx.==:l1a] and ARR[solvefor][:,mhidx.==:l1b]=ARR[:l1b][:,mhidx.==:l1b]
if solvefor :l1b, then allelem = [mhidx.==:l1b] and ARR[solvefor][:,mhidx.==:l1a]=ARR[:l1a][:,mhidx.==:l1a]
if solvefor 1, then allelem = [mhidx.==2; mhidx.==3]
if solvefor 2, then allelem = [mhidx.==2] and ARR[solvefor][:,mhidx.==3]=ARR[3][:,mhidx.==3]
if solvefor 3, then allelem = [mhidx.==3] and ARR[solvefor][:,mhidx.==2]=ARR[2][:,mhidx.==2]

# `activehypo` in example mh=[1.0;0.5;0.5]
sfidx=1, mhidx=2:  ah = [1;2]
sfidx=1, mhidx=3:  ah = [1;3]
sfidx=2, mhidx=2:  ah = [1;2]
sfidx=2, mhidx=3:  2 should take a value from 3 or nullhypo
sfidx=3, mhidx=2:  3 should take a value from 2 or nullhypo
sfidx=3, mhidx=3:  ah = [1;3]

# `activehypo` in example mh=[1.0;0.33;0.33;0.34]
sfidx=1, mhidx=2:  ah = [1;2]
sfidx=1, mhidx=3:  ah = [1;3]
sfidx=1, mhidx=4:  ah = [1;4]

sfidx=2, mhidx=2:  ah = [1;2]
sfidx=2, mhidx=3:  2 should take a value from 3 or nullhypo
sfidx=2, mhidx=4:  2 should take a value from 4 or nullhypo

sfidx=3, mhidx=2:  3 should take a value from 2 or nullhypo
sfidx=3, mhidx=3:  ah = [1;3]
sfidx=3, mhidx=4:  3 should take a value from 4 or nullhypo

sfidx=4, mhidx=2:  4 should take a value from 2 or nullhypo
sfidx=4, mhidx=3:  4 should take a value from 3 or nullhypo
sfidx=4, mhidx=4:  ah = [1;4]
```

Also keeping the default case documented:
```
# the default case where mh==nothing
# equivalent to mh=[1;1;1] # assuming 3 variables
sfidx=1, allelements=allidx[nhidx.==0], activehypo=(0,[1;])
sfidx=2, allelements=allidx[nhidx.==0], activehypo=(0,[2;])
sfidx=3, allelements=allidx[nhidx.==0], activehypo=(0,[3;])
```

TODO still need to compensate multihypo case for user nullhypo addition.
"""
function _prepareHypoRecipe!(
  mh::Categorical,
  maxlen::Int,
  sfidx::Int,
  lenXi::Int,
  isinit::Vector{Bool} = ones(Bool, lenXi),
  nullhypo::Real = 0,
)
  #
  allelements = Vector{Vector{Int}}()
  activehypo = Vector{Tuple{Int, Vector{Int}}}()
  mhidx = Vector{Int}()

  allidx = 1:maxlen
  allmhp, certainidx, uncertnidx = getHypothesesVectors(mh.p)

  # select only hypotheses that can be used (ie variables have been initialized)
  @assert !(sum(isinit) == 0 && sfidx == certainidx) # cannot init from nothing for any hypothesis

  mhh = if sum(isinit) < lenXi - 1
    @assert isLeastOneHypoAvailable(sfidx, certainidx, uncertnidx, isinit)
    @info "not all hypotheses initialized, but at least one available -- see #427"
    mhp = deepcopy(mh.p)
    suppressmask = isinit .== false
    suppressmask[sfidx] = false
    mhp[suppressmask] .= 0.0
    mhp ./= sum(mhp)
    Categorical(mhp)
  else
    mh
  end

  # FIXME consolidate with addEntropyOnManifolds approach in `computeAcrossHypothesis!`
  # prepend for the mhidx=0, bad-init-null-hypothesis case (if solving a fractional variable)
  mhh = if sfidx in uncertnidx
    nhw = (length(uncertnidx) + 1)
    nmhw = [1 / nhw; length(uncertnidx) / nhw * mhh.p]
    nmhw ./= sum(nmhw) # renormalize (should not be necessary)
    Categorical(nmhw)
  else
    mhh
  end

  # prep mm-nultihypothesis selection values
  mhidx = rand(mhh, maxlen)  # selection of which hypothesis is correct
  pidx = 0
  if sfidx in uncertnidx
    # shift down to get mhidx=0 case
    mhidx .-= 1
    pidx = -1
  end

  sfincer = sfidx in certainidx
  for pval in mhh.p
    pidx += 1
    pidxincer = pidx in certainidx # ??
    # permutation vectors for later computation
    iterarr = allidx[mhidx .== pidx]
    iterah = Int[]
    if !pidxincer && sfincer && pidx != 0 # 1e-15 <= pval && mh.p[sfidx] < 1e-10  # proxy for sfidx in certainidx
      # solve for one of the certain variables containing uncertain hypotheses in others
      iterah = sort(union(certainidx, pidx)) # sort([sfidx;pidx])
    # DONE -- supports n-ary factors in multihypo mode
    elseif (pidxincer && !sfincer || sfidx == pidx) && pidx != 0 # pval < 1e-15 && mh.p[sfidx] >= 1e-10
      # solve for one of the uncertain variables
      iterah = sort(union(certainidx, sfidx)) # sort([sfidx;pidx])
    # EXPERIMENTAL -- support more than binary factors in multihypo mode
    elseif pidxincer && sfincer && pidx != 0 # pval < 1e-15 && mh.p[sfidx] < 1e-10
      iterarr = Int[]
      iterah = Int[] # may be moot anyway, but double check first
    elseif !pidxincer && !sfincer && pidx != 0 # pval >= 1e-15 && mh.p[sfidx] >= 1e-10
      iterah = uncertnidx #allmhp[mh.p .> 1e-15]
    elseif pidx == 0
      # nullhypo only take values from self sfidx, might add entropy later (for bad init case)
      iterah = [sfidx;]
    else
      error("Unknown hypothesis case, got sfidx=$(sfidx) with mh.p=$(mh.p), pidx=$(pidx)")
    end
    push!(allelements, iterarr)
    push!(activehypo, (pidx, iterah))
  end

  # # retroactively add nullhypo case (the 0 case)
  # if sfidx in uncertnidx
  #   #
  # end

  return HypoRecipe(; certainidx, allelements, activehypo, mhidx)
  # hyporecipe::NamedTuple
  # return (; certainidx, allelements, activehypo, mhidx)
end

function _prepareHypoRecipe!(
  mh::Nothing,
  maxlen::Int,
  sfidx::Int,
  lenXi::Int,
  isinit::Vector{Bool} = ones(Bool, lenXi),
  nullhypo::Real = 0,
)
  #
  # FIXME, consolidate with the general multihypo case

  # the default case where mh==nothing
  # equivalent to mh=[1;1;1] # assuming 3 variables
  # sfidx=1, allelements=allidx[nhidx.==0], activehypo=(0,[1;])

  #
  allelements = Vector{Vector{Int}}()
  activehypo = Vector{Tuple{Int,Vector{Int}}}()

  # TODO add cases where nullhypo occurs, see DFG #536, and IIF #237
  nmhw = [nullhypo; (1 - nullhypo)]
  nhh = Categorical(nmhw)

  # prep mmultihypothesis selection values
  # mhidx = Int[]
  # NOTE, must do something special to get around Categorical([0;10]) error
  # selection of which hypothesis is correct
  mhidx = nullhypo == 0 ? ones(Int, maxlen) : (rand(nhh, maxlen) .- 1)

  allidx = 1:maxlen
  certainidx = 1:lenXi
  # zero is nullhypo case, 1 is first sfidx variable
  nullarr = allidx[mhidx .== 0]
  # mhidx == 1 case is regular -- this will be all elements if nullhypo=0.0
  reguarr = allidx[mhidx .!= 0]
  pidxAll = [0; certainidx]
  for pidx in pidxAll
    if pidx == 0
      # elements that occur during nullhypo active
      push!(allelements, nullarr)
      push!(activehypo, (pidx, [sfidx;]))
    elseif pidx == 1
      # elements that occur during regular hypothesis true
      push!(allelements, reguarr)
      push!(activehypo, (pidx, certainidx))
    else
      # all remaining collections are empty (part of multihypo support)
      push!(allelements, Int[])
      push!(activehypo, (pidx, Int[]))
    end
  end

  return HypoRecipe(; certainidx, allelements, activehypo, mhidx)
  # return hyporecipe::NamedTuple
  # return (; certainidx, allelements, activehypo, mhidx)
end

#
