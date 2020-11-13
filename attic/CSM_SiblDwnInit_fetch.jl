# Downward initialization through a cascading model
# these functions need to be consolidated refactored and deprecated accordingly
# see #910


## ============================================================================================
## Down init sibling priority order
## ============================================================================================


"""
$SIGNATURES

Test waiting order between siblings for cascading downward tree initialization.

Notes
- State machine function 8j.

DevNotes
- FIXME, this guy is never called during WIP of 459 dwnMsg consolidation
- FIXME, something wrong with CSM sequencing, https://github.com/JuliaRobotics/IncrementalInference.jl/issues/602#issuecomment-682114232
- This might be replaced with 4-stroke tree-init.
"""
function dwnInitSiblingWaitOrder_StateMachine(csmc::CliqStateMachineContainer)

  prnt_ = getParent(csmc.tree, csmc.cliq)
  if 0 == length(prnt_)
    # go to 7
    return determineCliqNeedDownMsg_StateMachine
  end
  prnt = prnt_[1]

  opt = getSolverParams(csmc.cliqSubFg) # csmc.dfg

  # now get the newly computed message from the appropriate container
  # make sure this is a pull model #674 (pull msg from source/prnt)
  # FIXME must be consolidated as part of #459
  # dwinmsgs = getfetchCliqueInitMsgDown(getCliqueData(prnt), from=:dwnInitSiblingWaitOrder_StateMachine)
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys_ = collect(keys(dwinmsgs.belief))
  infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, dwinmsgs keys=$(dwnkeys_)")

  # add downward belief prop msgs
  msgfcts = addMsgFactors!(csmc.cliqSubFg, dwinmsgs, DownwardPass)
  ## update solveable dims with newly avaiable init information
  # determine if more info is needed for partial
  sdims = getCliqVariableMoreInitDims(csmc.cliqSubFg, csmc.cliq)
  infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, sdims=$(sdims)")
  updateCliqSolvableDims!(csmc.cliq, sdims, csmc.logger)

  _dbgCSMSaveSubFG(csmc, "fg_DWNCMN_8j")

  ## FIXME, if this new, and all sibling clique's solvableDim are 0, then go back to waitChangeOnParentCondition_StateMachine

  # go to 8o.i
  return testDirectDwnInit_StateMachine
end

"""
    $SIGNATURES

Can this clique initialize directly from available down message info?

Notes
- State machine function nr. 8o.i
- Assume must have parent since this only occurs after 8j.
"""
function testDirectDwnInit_StateMachine(csmc::CliqStateMachineContainer)

  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys_ = collect(keys(dwinmsgs.belief))
  # NOTE, only use separators, not all parent variables (DF ???)
  # dwnkeys_ = lsf(csmc.cliqSubFg, tags=[:DOWNWARD_COMMON;]) .|> x->ls(csmc.cliqSubFg, x)[1]
  # @assert length(intersect(dwnkeys, dwnkeys_)) == length(dwnkeys) "split dwnkeys_ is not the same, $dwnkeys, and $dwnkeys_, separators: $(getCliqSeparatorVarIds(csmc.cliq))"

  # priorize solve order for mustinitdown with lowest dependency first
  # follow example from issue #344
  # mustwait = false
  if length(intersect(dwnkeys_, getCliqSeparatorVarIds(csmc.cliq))) == 0
    # can directly use DOWNWARD_COMMON
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, no can do, must wait for siblings to update parent first.")
    # mustwait = true

    # go to 8o.iv
    return sibsDwnPriorityInit_StateMachine
  end

  # go to 8o.ii
  return testDelayOrderDwnInit_StateMachine
end

"""
    $SIGNATURES

Return true if this clique's down init should be delayed on account of prioritization among sibling separators.

Notes
- State machine function nr. 8o.ii
- process described in issue #344

Dev Notes
- not priorizing order yet (TODO), just avoiding unsolvables at this time.
- Very closely related to 8o.iii -- refactor likely (NOTE).
- should precompute `allinters`.
- # FIXME ON FIRE, this is doing `getCliqueStatus` of siblings without following channels.
  - NOTE This is only used in old tree-init, not part of active regular solve code.
  - See #954
"""
function testDelayOrderDwnInit_StateMachine(csmc::CliqStateMachineContainer)

  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys = collect(keys(dwinmsgs.belief))
  
  tree = csmc.tree
  cliq = csmc.cliq
  logger = csmc.logger
  
  # when is a cliq upsolved
  solvedstats = Symbol[:upsolved; :marginalized; :uprecycled]

  # safety net double check
  cliqst = getCliqueStatus(cliq)
  if cliqst in solvedstats
    infocsm(csmc, "getSiblingsDelayOrder -- clique status should not be here with a solved cliqst=$cliqst")
    # go to 8o.iii
    return testPartialNeedsDwnInit_StateMachine
  end

  # get siblings separators
  sibs = getCliqSiblings(tree, cliq, true)
  ids = map(s->s.index, sibs)
  len = length(sibs)
  sibidx = collect(1:len)[ids .== cliq.index][1]
  seps = getCliqSeparatorVarIds.(sibs)
  lielbls = setdiff(ids, cliq.index)
  # get intersect matrix of siblings (should be exactly the same across siblings' csm)
  allinters = Array{Int,2}(undef, len, len)
  dwninters = Vector{Int}(undef, len)
  infocsm(csmc, "getSiblingsDelayOrder -- number siblings=$(len), sibidx=$sibidx")

  # sum matrix with all "up solved" rows and columns eliminated
  fill!(allinters, 0)
  for i in 1:len
    for j in i:len
      if i != j
        allinters[i,j] = length(intersect(seps[i],seps[j]))
      end
    end
    dwninters[i] = length(intersect(seps[i], dwnkeys))
  end

  # sum "across/over" rows, then columns (i.e. visa versa "along" columns then rows)
  rows = sum(allinters, dims=1)
  cols = sum(allinters, dims=2)

  infocsm(csmc, "getSiblingsDelayOrder -- allinters=$(allinters), getSiblingsDelayOrder -- rows=$(rows), getSiblingsDelayOrder -- rows=$(cols)")

  # is this clique a non-zero row -- i.e. sum across columns? if not, no further special care needed
  if cols[sibidx] == 0
    infocsm(csmc, "getSiblingsDelayOrder -- cols[sibidx=$(sibidx))] == 0, no special care needed")
    # go to 8o.iii
    return testPartialNeedsDwnInit_StateMachine
  end

  # now determine if initializing from below or needdownmsg
  if cliqst in Symbol[:needdownmsg;]
    # be super careful about delay (true) vs pass (false) at this point -- might be partial too TODO
    # return true if delay beneficial to initialization accuracy

    # find which siblings this cliq epends on
    symm = allinters + allinters'
    maskcol = 0 .< symm[:,sibidx]
    # lenm = length(maskcol)
    stat = Vector{Symbol}(undef, len)
    stillbusymask = fill(false, len)

    # get each sibling status (entering atomic computation segment -- until wait command)
    stat .= getCliqueStatus.(sibs) #[maskcol]

    ## (long down chain case)
    # need different behaviour when all remaining siblings are blocking with :needdownmsg
    remainingmask = stat .== :needdownmsg
    if sum(remainingmask) == length(stat)
      infocsm(csmc, "getSiblingsDelayOrder -- all blocking: sum(remainingmask) == length(stat), stat=$stat")
      # pick sibling with most overlap in down msgs from parent
      # list of similar length siblings
      candidates = dwninters .== maximum(dwninters)
      if candidates[sibidx]
        # must also pick minimized intersect with other remaing siblings
        maxcan = collect(1:len)[candidates]
        infocsm(csmc, "getSiblingsDelayOrder -- candidates=$candidates, maxcan=$maxcan, rows=$rows")
        if rows[sibidx] == minimum(rows[maxcan])
          infocsm(csmc, "getSiblingsDelayOrder -- FORCE DOWN INIT SOLVE ON THIS CLIQUE: $(cliq.index), $(getLabel(cliq))")
          # go to 8o.iii
          return testPartialNeedsDwnInit_StateMachine
        end
      end
      infocsm(csmc, "getSiblingsDelayOrder -- not a max and should block")
      # go to 8o.iv
      return sibsDwnPriorityInit_StateMachine
    end

    # still busy solving on branches, so potential to delay
    for i in 1:len
      stillbusymask[i] = maskcol[i] && !(stat[i] in solvedstats)
    end
    infocsm(csmc, "getSiblingsDelayOrder -- busy solving: maskcol=$maskcol, stillbusy=$stillbusymask")

    # Too blunt -- should already have returned false by this point perhaps
    if 0 < sum(stillbusymask)
      # yes something to delay about
      infocsm(csmc, "getSiblingsDelayOrder -- yes delay, stat=$stat symm=$symm")

      # go to 8o.iv
      return sibsDwnPriorityInit_StateMachine
    end
  end

  infocsm(csmc, "getSiblingsDelayOrder -- default will not delay")
  # carry over default from partial init process

  # go to 8o.iii
  return testPartialNeedsDwnInit_StateMachine
end



"""
    $SIGNATURES


Return true if both, i.) this clique requires more downward information, ii.) more
downward message information could potentially become available.

Notes
- State machine function nr. 8o.iii
- Delay initialization to the last possible moment.

Dev Notes:
- # FIXME ON FIRE, this is doing `getCliqueStatus` of siblings without following channels.
  - NOTE This is only used in old tree-init, not part of active regular solve code.
  - See #954
- Determine clique truely isn't able to proceed any further:
  - should be as self reliant as possible (using clique's status as indicator)
  - OBSOLETE ?? 
    - change status to :mustinitdown if have only partial beliefs so far:
    - combination of status, while partials belief siblings are not :mustinitdown
"""
function testPartialNeedsDwnInit_StateMachine(csmc::CliqStateMachineContainer)
  #
  prnt = getParent(csmc.tree, csmc.cliq)[1]
  dwinmsgs = fetchDwnMsgConsolidated(prnt)

  tree = csmc.tree
  cliq = csmc.cliq
  logger = csmc.logger

  # which incoming messages are partials
  hasPartials = Dict{Symbol, Int}()
  for (sym, tmsg) in dwinmsgs.belief
    # assuming any down message per label that is not partial voids further partial consideration
    if sum(tmsg.inferdim) > 0
      if !haskey(hasPartials, sym)
        hasPartials[sym] = 0
      end
      hasPartials[sym] += 1
    end
  end
  partialKeys = collect(keys(hasPartials))

  ## determine who might be able to help init this cliq
  # check sibling separator sets against this clique's separator
  sibs = getCliqSiblings(tree, cliq)

  infocsm(csmc, "getCliqSiblingsPartialNeeds -- CHECK PARTIAL")
  # identify which cliques might have useful information
  localsep = getCliqSeparatorVarIds(cliq)
  seps = Dict{Int, Vector{Symbol}}()
  for si in sibs
    # @show getLabel(si)
    mighthave = intersect(getCliqSeparatorVarIds(si), localsep)
    if length(mighthave) > 0
      seps[si.index] = mighthave
      if getCliqueStatus(si) in [:initialized; :null; :needdownmsg]
        # partials treated special -- this is slightly hacky
        if length(intersect(localsep, partialKeys)) > 0 && length(mighthave) > 0
          # this sibling might have info to delay about
          setCliqueDrawColor!(cliq,"magenta")
          # go to 8o.iv
          return sibsDwnPriorityInit_StateMachine
        end
      end
    end
  end
  # determine if those cliques will / or will not be able to provide more info
  # when does clique change to :mustinitdown
  # default

  # go to 8e.ii.
  return tryDwnInitCliq_StateMachine
end

    
    
"""
$SIGNATURES
    
    
Return true there is no other sibling that will make progress.

Notes
- State machine function nr. 8o.iv
- Relies on sibling priority order with only one "currently best" option that will force progress in global upward inference.
- Return false if one of the siblings is still busy

DevNotes
- # FIXME ON FIRE, calling getCliqueStatus of siblings directly
  - Only used for older tree init, not part of regular up-down solving code
  - see #954
- Best is likely if parent actually made determination on which child to solve first in this case.
"""
function sibsDwnPriorityInit_StateMachine(csmc::CliqStateMachineContainer)

  tree = csmc.tree
  cliq = csmc.cliq
  prnt_ = getParent(tree, cliq)
  prnt = prnt_[1]

  dwinmsgs = fetchDwnMsgConsolidated(prnt)
  dwnkeys_ = collect(keys(dwinmsgs.belief))

  solvord = getCliqSiblingsPriorityInitOrder( csmc.tree, prnt, dwinmsgs, csmc.logger )

  # noOneElse = areSiblingsRemaingNeedDownOnly(csmc.tree, csmc.cliq)
  noOneElse = true 
  stillbusylist = [:null; :initialized;]
  if 0 < length(prnt_)
    for si in getChildren(tree, prnt)
      # are any of the other siblings still busy?
      if si.index != cliq.index && getCliqueStatus(si) in stillbusylist
        noOneElse = false
      end
    end
  end
  # nope, everybody is waiting for something to change -- proceed with forcing a cliq solve

  infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, $(prnt.index), $noOneElse, solvord = $solvord")


  if csmc.cliq.index != solvord[1]
    # TODO, is this needed? fails hex if included with correction :needdowninit --> :needdownmsg
    # if dwinmsgs.status == :initialized && getCliqueStatus(csmc.cliq) == :needdownmsg  # used to be :needdowninit
    #   # go to 7e
    #   return slowWhileInit_StateMachine
    # end

    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, must wait on change.")
    # remove all message factors
    # remove msg factors previously added
    fctstorm = deleteMsgFactors!(csmc.cliqSubFg, [:DOWNWARD_COMMON])
    infocsm(csmc, "8j, dwnInitSiblingWaitOrder_StateMachine, removing factors $fctstorm")

    # go to 8c
    return waitChangeOnParentCondition_StateMachine
  end

  # go to 8e.ii.
  return tryDwnInitCliq_StateMachine
end








