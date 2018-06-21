# Second iteration of explicitly exloring the marginalization of discrete variables onto the continuous state space.
# Although this code is still excessive and messy, this is a significant feature expansion; and future versions will
# generalize the marginalization process to allow for implicit hypothesis exploration.  The messy explicit version of code is
# intended to help develop the required generalistic unit tests.  The unit tests will then be validated, frozen and uesd to
# confirm future "algebraic" marginalization (implicit) versions operate correctly.
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

# `activehypo` in example mh=[0;0.5;0.5]
sfidx=1, mhidx=2:  ah = [1;2]
sfidx=1, mhidx=3:  ah = [1;3]
sfidx=2, mhidx=2:  ah = [1;2]
sfidx=2, mhidx=3:  2 should take a value from 3
sfidx=3, mhidx=2:  3 should take a value from 2
sfidx=3, mhidx=3:  ah = [1;3]

# `activehypo` in example mh=[0;0.33;0.33;0.34]
sfidx=1, mhidx=2:  ah = [1;2]
sfidx=1, mhidx=3:  ah = [1;3]
sfidx=1, mhidx=4:  ah = [1;4]
...
sfidx=2, mhidx=3:  2 should take a value from 3
```
"""
function assembleHypothesesElements!(
            # allelements::Array,
            # activehypo::Array,
            mh::Categorical,
            maxlen::Int,
            sfidx::Int,
            # mhidx::Vector{<:Int},
            lenXi::Int  )
  #
  allelements = []
  activehypo = []
  mhidx = Int[]

  mh.p
  allidx = 1:maxlen
  allmhp = 1:length(mh.p)
  certainidx = allmhp[mh.p .< 1e-10]

  if mh != nothing
    # If present, prep mmultihypothesis selection values
    mhidx = rand(mh, maxlen) # selection of which hypothesis is correct
  end

  # this is not going to work? sfidx could be anything
  if mh.p[sfidx] < 1e-10
    pidx = 0
    for pval in mh.p
      pidx += 1
      if pval > 1e-10
        iterarr = allidx[mhidx .== pidx]
        push!(allelements, iterarr)
        iterah = sort([sfidx;pidx]) # TODO -- currently only support binary factors in multihypo mode
        push!(activehypo, (pidx, iterah))
      end
    end
  elseif mh.p[sfidx] >= 1e-10
    pidx = 0
    for pval in mh.p
      pidx += 1
      # must still include cases where sfidx != pidx
      ## TODO -- Maybe a mistake with iterah variables in these cases?
      if pval < 1e-10
        iterarr = allidx[mhidx .== pidx]
        push!(allelements, iterarr)
        iterah = sort([sfidx;pidx]) # TODO -- currently only support binary factors in multihypo mode
        push!(activehypo, (pidx, iterah))
      elseif pval > 1e-10 && sfidx == pidx
        iterarr = allidx[mhidx .== pidx]
        push!(allelements, iterarr)
        iterah = allmhp[mh.p .> 1e-10]
        # @show iterah = sort(union([pidx;], allmhp[mh.p .< 1e-10]))
        push!(activehypo, (pidx,iterah))
      elseif pval > 1e-10 && sfidx != pidx
        iterarr = allidx[mhidx .== pidx]
        push!(allelements, iterarr)
        iterah = allmhp[mh.p .> 1e-10]
        push!(activehypo, (pidx,iterah))
      else
        error("assembleHypothesesElements! mh.p[sfidx=$(sfidx)] >= 1e-10 is missing a case: pval=$pval")
      end
    end
  else
    error("Unknown hypothesis case, got sfidx=$(sfidx) with mh.p=$(mh.p)")
  end

  return certainidx, allelements, activehypo, mhidx
end
function assembleHypothesesElements!(
            # allelements::Array,
            # activehypo::Array,
            mh::Void,
            maxlen::Int,
            sfidx::Int,
            # mhidx,
            lenXi::Int  )
  #
  allelements = []
  activehypo = []
  mhidx = Int[]

  # error("assembleHypothesesElements!(..) -- Error in code design, refactor of general multihypothesis situations required if you arrived here.")
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
