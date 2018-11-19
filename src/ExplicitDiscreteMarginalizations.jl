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
    $(SIGNATURES)

This function explicitly codes that marginalization of a discrete categorical selection variable for ambiguous data association situations.  Improved implementations should implicitly induce the same behaviour through summation (integration) when marginalizing any number of discrete variables.  This function populates `allelements` with particle indices associated with particular multihypothesis selection while `activehypo` simultaneously contains the hypothesis index and factor graph variables associated with that hypothesis selection.  The return value `certainidx` are the hypotheses that are not in question.

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
sfidx=2, mhidx=3:  2 should take a value from 3
sfidx=3, mhidx=2:  3 should take a value from 2
sfidx=3, mhidx=3:  ah = [1;3]

# `activehypo` in example mh=[1.0;0.33;0.33;0.34]
sfidx=1, mhidx=2:  ah = [1;2]
sfidx=1, mhidx=3:  ah = [1;3]
sfidx=1, mhidx=4:  ah = [1;4]

sfidx=2, mhidx=2:  ah = [1;2]
sfidx=2, mhidx=3:  2 should take a value from 3
sfidx=2, mhidx=4:  2 should take a value from 4

sfidx=3, mhidx=2:  3 should take a value from 2
sfidx=3, mhidx=3:  ah = [1;3]
sfidx=3, mhidx=4:  3 should take a value from 4

sfidx=4, mhidx=2:  4 should take a value from 2
sfidx=4, mhidx=3:  4 should take a value from 3
sfidx=4, mhidx=4:  ah = [1;4]
```
"""
function assembleHypothesesElements!(
            mh::Categorical,
            maxlen::Int,
            sfidx::Int,
            lenXi::Int  )
  #
  allelements = []
  activehypo = []
  mhidx = Int[]

  allidx = 1:maxlen
  allmhp = 1:length(mh.p)
  certainidx = allmhp[mh.p .== 0.0]  # TODO remove after gwp removed
  uncertnidx = allmhp[0.0 .< mh.p]

  # prep mmultihypothesis selection values
  mhidx = rand(mh, maxlen) # selection of which hypothesis is correct

  pidx = 0
  sfincer = sfidx in certainidx
  for pval in mh.p
    pidx += 1
    pidxincer = pidx in certainidx # ??
    # permutation vectors for later computation
    iterarr = allidx[mhidx .== pidx]
    iterah = Int[]
    if !pidxincer && sfincer # 1e-15 <= pval && mh.p[sfidx] < 1e-10  # proxy for sfidx in certainidx
      # solve for one of the certain variables containing uncertain hypotheses in others
      iterah = sort(union(certainidx, pidx)) # sort([sfidx;pidx])
      # DONE -- supports n-ary factors in multihypo mode
    elseif pidxincer && !sfincer || sfidx == pidx # pval < 1e-15 && mh.p[sfidx] >= 1e-10
      # solve for one of the uncertain variables
      iterah = sort(union(certainidx, sfidx)) # sort([sfidx;pidx])
      # EXPERIMENTAL -- support more than binary factors in multihypo mode
    elseif pidxincer && sfincer # pval < 1e-15 && mh.p[sfidx] < 1e-10
      iterarr = Int[]
      iterah = Int[] # may be moot anyway, but double check first
    elseif !pidxincer && !sfincer # pval >= 1e-15 && mh.p[sfidx] >= 1e-10
      iterah = uncertnidx #allmhp[mh.p .> 1e-15]
    else
      error("Unknown hypothesis case, got sfidx=$(sfidx) with mh.p=$(mh.p), pidx=$(pidx)")
    end
    push!(allelements, iterarr)
    push!(activehypo, (pidx,iterah))
  end

  # if mhidx == sfidx
  # ah = sort(union([sfidx;], certainidx))

  return certainidx, allelements, activehypo, mhidx
end
function assembleHypothesesElements!(
            mh::Nothing,
            maxlen::Int,
            sfidx::Int,
            lenXi::Int  )
  #
  allelements = []
  activehypo = []
  mhidx = Int[]

  allidx = 1:maxlen
  certainidx = 1:lenXi
  doneall = false
  for i in certainidx
    if !doneall
      push!(allelements, allidx)
      push!(activehypo, (i,certainidx))
      doneall = true
    else
      push!(allelements, Int[])
      push!(activehypo, (i,Int[]))
    end
  end
  return certainidx, allelements, activehypo, mhidx # certainidx = allhp
end







#
