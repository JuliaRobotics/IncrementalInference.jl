
"""
    $SIGNATURES

Based on a push model from child cliques that should have already completed their computation.
"""
function getCliqInitUpMsgs(cliq::Graphs.ExVertex)
  getData(cliq).upInitMsgs
end

function setCliqUpInitMsgs!(cliq::Graphs.ExVertex, childid::Int, msg::Dict)
  getData(cliq).upInitMsgs[childid] = msg
end

function isCliqInitialized(cliq::Graphs.ExVertex)::Bool
  return getData(cliq).initialized in [:initialized; :upsolved]
end

function isCliqUpSolved(cliq::Graphs.ExVertex)::Bool
  return getData(cliq).initialized == :upsolved
end


"""
    $SIGNATURES

Return the most likely  ordering for initializing factor (assuming up solve
sequence).

Notes:
- sorts id for increasing number of connected factors.
"""
function getCliqInitVarOrderUp(cliq::Graphs.ExVertex)
  # rules to explore dimension from one to the other?

  # get all variable ids and number of associated factors
  allids = getCliqAllVarIds(cliq)
  nfcts = getCliqNumAssocFactorsPerVar(cliq)

  # get priors and singleton message variables (without partials)
  prids = getCliqVarIdsPriors(cliq, getCliqAllVarIds(cliq), false)

  # get current up msgs in the init process (now have all singletons)
  upmsgs = getCliqInitUpMsgs(cliq)
  upmsgids = collect(keys(upmsgs))

  # all singleton variables
  singids = union(prids, upmsgids)

  # add msg marginal prior (singletons) to number of factors
  for msid in upmsgids
    nfcts[msid .== allids] .+= 1
  end

  # sort permutation order for increasing number of factor association
  nfctsp = sortperm(nfcts)
  sortedids = allids[nfctsp]

  # organize the prior variables separately with asceding factor count
  initorder = zeros(Int, 0)
  for id in sortedids
    if id in singids
      push!(initorder, id)
    end
  end
  # in ascending order of number of factors
  for id in sortedids
    if !(id in initorder)
      push!(initorder, id)
    end
  end
  return initorder
end

function notifyCliqUpInitStatus!(cliq::Graphs.ExVertex, status::Symbol)
  cd = getData(cliq)
  cd.initialized = status
  if isready(cd.initUpChannel)
    tkst = take!(cd.initUpChannel)
    @info "dumping stale cliq=$(cliq.index) status message $(tkst), replacing with $(status)"
  end
  put!(cd.initUpChannel, status)
end

function notifyCliqDownInitStatus!(cliq::Graphs.ExVertex, status::Symbol)
  @info "$(current_task()) Clique $(cliq.index), notify down init status=$(status)"
  cd = getData(cliq)
  cd.initialized = status
  if isready(cd.initDownChannel)
    @info "dumping stale cliq=$(cliq.index) status message $(take!(cd.initDownChannel)), replacing with $(status)"
  end
  put!(cd.initDownChannel, status)
end

"""
    $SIGNATURES

Return true if clique has completed the local upward direction inference procedure.
"""
isUpInferenceComplete(cliq::Graphs.ExVertex) = getData(cliq).upsolved

function areCliqVariablesAllInitialized(fgl::FactorGraph, cliq::Graphs.ExVertex)
  allids = getCliqAllVarIds(cliq)
  isallinit = true
  for vid in allids
    var = getVert(fgl, vid, api=localapi)
    isallinit &= isInitialized(var)
  end
  isallinit
end

"""
   $SIGNATURES

Determine if this `cliq` has been fully initialized and child cliques have completed their full upward inference.
"""
function isCliqReadyInferenceUp(fgl::FactorGraph, tree::BayesTree, cliq::Graphs.ExVertex)
  isallinit = areCliqVariablesAllInitialized(fgl, cliq)

  # check that all child cliques have also completed full up inference.
  for chl in getChildren(tree, cliq)
    isallinit &= isUpInferenceComplete(chl)
  end
  return isallinit
end

"""
    $SIGNATURES

Blocking call until `cliq` upInit processes has arrived at a result.
"""
function getCliqInitUpResultFromChannel(cliq::Graphs.ExVertex)
  status = take!(getData(cliq).initUpChannel)
  @info "$(current_task()) Clique $(cliq.index), dumping initUpChannel status, $status"
  return status
end

"""
    $SIGNATURES

Return `::Symbol` status a particular clique is in, with specific regard to solution
or numerical initialization status:
- :needdownmsg
- :upsolved
- :downsolved
- :initialized
- :marginalized
- :null

Notes:
- `:null` represents the first uninitialized state of a cliq.
"""
getCliqStatusUp(cliq::Graphs.ExVertex)::Symbol = getData(cliq).initialized
getCliqStatus(cliq::Graphs.ExVertex)::Symbol = getCliqStatusUp(cliq)

"""
    $SIGNATURES

Set up initialization or solve status of this `cliq`.
"""
function setCliqStatus!(cliq::Graphs.ExVertex, status::Symbol)
  getData(cliq).initialized = status
end




"""
    $SIGNATURES

Return true if all variables in clique are considered marginalized (and initialized).
"""
function areCliqVariablesAllMarginalized(subfg::FactorGraph, cliq::Graphs.ExVertex)
  for vid in getCliqAllVarIds(cliq)
    vert = getVert(subfg, vid)
    if !isMarginalized(vert) || !isInitialized(vert)
      return false
    end
  end
  return true
end


"""
    $SIGNATURES

Set all Bayes (Junction) tree cliques that have all marginalized and initialized variables.
"""
function setTreeCliquesMarginalized!(fgl::FactorGraph, tree::BayesTree)
  for (cliid, cliq) in tree.cliques
    if areCliqVariablesAllMarginalized(fgl, cliq)
      # need to set the upward messages
      msgs = prepCliqInitMsgsUp(fgl, cliq)
      setUpMsg!(cliq, msgs)

      prnt = getParent(tree, cliq)
      if length(prnt) > 0
        # THIS IS FOR INIT PASSES ONLY
        setCliqUpInitMsgs!(prnt[1], cliq.index, msgs)
      end

      setCliqStatus!(cliq, :marginalized)
      setCliqDrawColor(cliq, "blue")
    end
  end
  nothing
end





"""
    $SIGNATURES

Block the thread until child cliques of `prnt::Graphs.ExVertex` have finished
attempting upward initialization -- i.e. have status result.
Return `::Dict{Symbol}` indicating whether next action that should be taken
for each child clique.

Notes:
- See status options at `getCliqStatusUp(..)`.
- Can be called multiple times
"""
function blockCliqUntilChildrenHaveUpStatus(tree::BayesTree,
                                            prnt::Graphs.ExVertex  )::Dict{Int, Symbol}
  #
  ret = Dict{Int, Symbol}()
  chlr = getChildren(tree, prnt)
  # @sync begin
    for ch in chlr
      # either wait to fetch new result, or report or result
      chst = getCliqStatusUp(ch)
      @info "$(current_task()) Clique $(prnt.index), child $(ch.index) status is $(chst), isready(initUpCh)=$(isready(getData(ch).initUpChannel))."
      # if chst != :null && !isready(getData(ch).initUpChannel)
      #   ret[ch.index] = chst
      # else
        ret[ch.index] = fetch(getData(ch).initUpChannel)
      # end
    end
  # end
  return ret
end

"""
    $SIGNATURES

Wait here if all siblings and the parent status are `:needdownmsg`.
Return true when parent is `:initialized` after all were `:needdownmsg`

Notes
- used for regulating long need down message chains.
- exit strategy is parent becomes status `:initialized`.
"""
function blockCliqSiblingsParentNeedDown(tree::BayesTree,
                                         cliq::Graphs.ExVertex)
  #
  # ret = Dict{Int, Symbol}()
  prnt = getParent(tree, cliq)
  allneeddwn = true
  if length(prnt) > 0
    prstat = getCliqStatus(prnt[1])
    if prstat == :needdownmsg
      for ch in getChildren(tree, prnt[1])
        chst = getCliqStatusUp(ch)
        if chst != :needdownmsg
          allneeddwn = false
        end
      end
      if allneeddwn
        @warn "$(current_task()) Clique $(cliq.index), block since all siblings/parent needdownmsg."
        prtmsg = fetch(getData(prnt[1]).initDownChannel)
        @info "$(current_task()) Clique $(prnt[1].index), blockCliqSiblingsParentNeedDown -- after fetch $prstat, $prtmsg"
        if prtmsg == :initialized
          return true
        else
          @warn "$(current_task()) Clique $(prnt[1].index), maybe clear down init message $prtmsg"
          # take!(getData(prnt[1]).initDownChannel)
        end
      end
    end
  end
  return false
end

"""
    $SIGNATURES

Cycle through var order and initialize variables as possible in `subfg::FactorGraph`.
Return true if something was updated.

Notes:
- assumed `subfg` is a subgraph containing only the factors that can be used.
  - including the required up or down messages
- intended for both up and down initialization operations.
"""
function cycleInitByVarOrder!(subfg::FactorGraph, varorder::Vector{Int})::Bool
  @info "cycleInitByVarOrder! -- varorder=$(varorder)"
  retval = false
  count = 1
  while count > 0
    count = 0
    for vid in varorder
      var = getVert(subfg, vid, api=localapi)
      isinit = isInitialized(var)
      @info "going for doautoinit!, var.label=$(var.label), isinit=$(isinit)"
      # TODO -- must use factors and values in cliq (assume subgraph?)
      doautoinit!(subfg, ExVertex[var;], api=localapi)
      if isinit != isInitialized(var)
        count += 1
        retval = true
      end
    end
  end
  @info "cycleInitByVarOrder!, retval=$(retval)"
  return retval
end

"""
    $SIGNATURES

Update `subfg::FactorGraph` according to internal computations for a full upsolve.
"""
function doCliqUpSolve!(subfg::FactorGraph,
                        tree::BayesTree,
                        cliq::Graphs.ExVertex  )::Symbol
  #
  csym = Symbol(getVert(subfg, getCliqFrontalVarIds(cliq)[1], api=localapi).label)
  approxCliqMarginalUp!(subfg, tree, csym, false)
  getData(cliq).upsolved = true
  return :upsolved
end

"""
    $SIGNATURES

Prepare the upward inference messages from clique to parent and return as `Dict{Symbol}`.
"""
function prepCliqInitMsgsUp(subfg::FactorGraph,
                             cliq::Graphs.ExVertex)::Dict{Symbol, BallTreeDensity}
  #
  # construct init's up msg to place in parent from initialized separator variables
  msg = Dict{Symbol, BallTreeDensity}()
  for vid in getCliqSeparatorVarIds(cliq)
    var = getVert(subfg, vid, api=localapi)
    if isInitialized(var)
      msg[Symbol(var.label)] = getKDE(var)
    end
  end
  return msg
end

function prepCliqInitMsgsUp(subfg::FactorGraph, tree::BayesTree, cliq::Graphs.ExVertex)::Dict{Symbol, BallTreeDensity}
  @warn "deprecated, use prepCliqInitMsgsUp(subfg::FactorGraph, cliq::Graphs.ExVertex) instead"
  prepCliqInitMsgsUp(subfg, cliq)
end

"""
    $SIGNATURES

Perform cliq initalization calculation based on current state of the tree and factor graph,
using upward message passing logic.

> NOTE WORK IN PROGRESS

Notes
- Return either of (:initialized, :upsolved, :needdownmsg, :badinit)
- must use factors in cliq only, ensured by using subgraph -- TODO general case.
"""
function doCliqAutoInitUp!(subfg::FactorGraph,
                           tree::BayesTree,
                           cliq::Graphs.ExVertex;
                           up_solve_if_able::Bool=true, )::Symbol
  #
  # init up msg has special procedure for incomplete messages
  status = :needdownmsg
  varorder = Int[]

  # get incoming clique up messages
  upmsgs = getCliqInitUpMsgs(cliq)

  # add incoming up messages as priors to subfg
  @info "$(current_task()) Clique $(cliq.index), doCliqAutoInitUp! -- adding up message factors"
  msgfcts = addMsgFactors!(subfg, upmsgs)

    # TEMP
    # writeGraphPdf(subfg, show=true)

  # attempt initialize if necessary
  if !areCliqVariablesAllInitialized(subfg, cliq)
    # structure for all up message densities computed during this initialization procedure.
    varorder = getCliqInitVarOrderUp(cliq)
    # do physical inits, ignore cycle return value
    @info "$(current_task()) Clique $(cliq.index), doCliqAutoInitUp! -- going for up cycle order"
    cycleInitByVarOrder!(subfg, varorder)
    @info "$(current_task()) Clique $(cliq.index), doCliqAutoInitUp! -- finished with up cycle order"
  end

  # check if all cliq vars have been initialized so that full inference can occur on clique
  if areCliqVariablesAllInitialized(subfg, cliq)
    @info "$(current_task()) Clique $(cliq.index), going for doCliqUpSolve!, clique status = $(status)"
    status = doCliqUpSolve!(subfg, tree, cliq)
  else
    @info "$(current_task()) Clique $(cliq.index), all varialbes not initialized, status = $(status)"
  end

  # construct init's up msg to place in parent from initialized separator variables
  @info "$(current_task()) Clique $(cliq.index), going to prepCliqInitMsgsUp"
  msg = prepCliqInitMsgsUp(subfg, cliq) # , tree

  # put the init result in the parent cliq.
  prnt = getParent(tree, cliq)
  @info "$(current_task()) Clique $(cliq.index), prnt = getParent(tree, cliq) = $(prnt)"
  if length(prnt) > 0
    # not a root clique
    @info "$(current_task()) Clique $(cliq.index), doCliqAutoInitUp! -- putting upinitmsg in prnt=$(prnt[1].index), with msgs for $(collect(keys(msg)))"
    setCliqUpInitMsgs!(prnt[1], cliq.index, msg)
  end

  # remove msg factors that were added to the subfg
  @info "$(current_task()) Clique $(cliq.index), doCliqAutoInitUp! -- removing up message factors, length=$(length(msgfcts))"
  deleteMsgFactors!(subfg, msgfcts)

  # set flags in clique for multicore sequencing
  @info "$(current_task()) Clique $(cliq.index), doCliqAutoInitUp! -- sending notification of up init status"
  notifyCliqUpInitStatus!(cliq, status)
  return status
end


"""
    $SIGNATURES

Initialization downward message passing is different from regular inference since
it is possible that none of the child cliq variables have been initialized.

Notes
- init msgs from child upward passes are individually stored in this `cliq`.
- fresh product of overlapping beliefs are calculated on each function call.
"""
function prepCliqInitMsgsDown!(fgl::FactorGraph,
                               tree::BayesTree,
                               cliq::Graphs.ExVertex )
  #
  @info "$(current_task()) Clique $(cliq.index), prepCliqInitMsgsDown!"
  # get the current messages stored in the parent
  currmsgs = getCliqInitUpMsgs(cliq)
  @info "$(current_task()) Clique $(cliq.index), msg keys=$(collect(keys(currmsgs)))"

  # check if any msgs should be multiplied together for the same variable
  msgspervar = Dict{Symbol, Vector{BallTreeDensity}}()
  for (cliqid, msgs) in currmsgs
    @show cliqid, length(msgs)
    for (msgsym, msg) in msgs
      if !haskey(msgspervar, msgsym)
        msgspervar[msgsym] = Vector{BallTreeDensity}()
      end
      push!(msgspervar[msgsym], msg)
    end
  end

  @info "$(current_task()) Clique $(cliq.index), keys with msgs=$(collect(keys(msgspervar)))"

  # use default allocated dict
  products = getData(cliq).downInitMsg
  # multiply multiple messages together
  for (msgsym, msgs) in msgspervar
    # check if this particular down message requires msgsym
    if haskey(fgl.IDs, msgsym)
      if length(msgspervar[msgsym]) > 1
        products[msgsym] = manifoldProduct(msgs, getManifolds(fgl, msgsym))
      else
        products[msgsym] = msgs[1]
      end
    else
      # not required, therefore remove from message to avoid confusion
      if haskey(products, msgsym)
        delete!(products, msgsym)
      end
    end
  end

  @info "$(current_task()) Clique $(cliq.index), product keys=$(collect(keys(products)))"
  return products
end

"""

Special function to do initialization in downward direction, assuming that not all
variables can be initialized.  Relies on outside down messages.

Notes:
- assumed this `cliq` is being initialized from a previous `:needdownmsg` status.
- will use all possible local factors of cliquq in initilization process
- similar to upward initialization, but uses different message structure
  - first draft assumes upward messages will not be used,
  - full up solve still required which explicitly depends on upward messages.
- TODO replace with nested 'minimum degree' type variable ordering
"""
function getCliqInitVarOrderDown(fgl::FactorGraph,
                                 cliq::Graphs.ExVertex,
                                 downmsgs::Dict{Symbol, BallTreeDensity}  )
  #
  allids = getCliqAllVarIds(cliq)
  # convert input downmsg var symbols to integers (also assumed as prior beliefs)
  # make sure ids are in the clique set, since parent may have more variables.
  dwnmsgsym = intersect(collect(keys(downmsgs)), collect(keys(fgl.IDs)))
  dwnmsgids =  map(x -> fgl.IDs[x], dwnmsgsym )
  dwnvarids = intersect(allids, dwnmsgids)

  # find any other prior factors (might have partials)
  prvarids = getCliqVarIdsPriors(cliq, allids, true)
  hassinglids = union(dwnvarids, prvarids)

  # Get all other variable factor counts
  nfcts = getCliqNumAssocFactorsPerVar(cliq)
  # add msg marginal prior (singletons) to number of factors
  for msid in dwnmsgids
    nfcts[msid .== allids] .+= 1
  end

  # sort permutation order for increasing number of factor association
  nfctsp = sortperm(nfcts)
  sortedids = allids[nfctsp]

  # all singleton variables
  singids = union(prvarids, dwnvarids)

  # organize the prior variables separately with asceding factor count
  initorder = zeros(Int, 0)
  for id in sortedids
    if id in singids
      push!(initorder, id)
    end
  end
  # sort remaining variables for increasing associated factors
  for id in sortedids
    if !(id in initorder)
      push!(initorder, id)
    end
  end

  # return variable order
  return initorder
end

"""
    $SIGNATURES

Modify the subgraph`::FactorGraph` to include `msgs` as priors that are used
during clique inference.

Notes
- May be used initialization or inference, in both upward and downward directions.

Related

`deleteMsgFactors!`
"""
function addMsgFactors!(subfg::FactorGraph,
                        msgs::Dict{Symbol, BallTreeDensity})::Vector{ExVertex}
  # add messages as priors to this sub factor graph
  msgfcts = Graphs.ExVertex[]
  svars = union(ls(subfg)...)
  mvid = getMaxVertId(subfg)
  for (msym, dm) in msgs
    if msym in svars
      # @show "adding down msg $msym"
      mvid += 1
      fc = addFactor!(subfg, [msym], Prior(dm), autoinit=false, uid=mvid)
      push!(msgfcts, fc)
    end
  end
  return msgfcts
end

function addMsgFactors!(subfg::FactorGraph,
                        msgs::Dict{Symbol, Vector{BallTreeDensity}})::Vector{ExVertex}
  # add messages as priors to this sub factor graph
  msgfcts = Graphs.ExVertex[]
  svars = union(ls(subfg)...)
  mvid = getMaxVertId(subfg)
  for (msym, dms) in msgs
    for dm in dms
      if msym in svars
        # @show "adding down msg $msym"
        mvid += 1
        fc = addFactor!(subfg, [msym], Prior(dm), autoinit=false, uid=mvid)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end

function addMsgFactors!(subfg::FactorGraph,
                        allmsgs::Dict{Int,Dict{Symbol, BallTreeDensity}})::Vector{Graphs.ExVertex}
  #
  allfcts = Graphs.ExVertex[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs)
    union!( allfcts, newfcts )
  end
  return allfcts
end

"""
    $SIGNATURES

Delete from the subgraph`::FactorGraph` prior belief `msgs` that could/would be used
during clique inference.

Related

`addMsgFactors!`
"""
function deleteMsgFactors!(subfg::FactorGraph, fcts::Vector{Graphs.ExVertex})
  for fc in fcts
    deleteFactor!(subfg, Symbol(fc.label))
  end
end


"""
    $SIGNATURES

Return true or false depending on whether child cliques are all up solved.
"""
function areCliqChildrenAllUpSolved(treel::BayesTree, prnt::Graphs.ExVertex)::Bool
  for ch in getChildren(treel, prnt)
    if !isCliqUpSolved(ch)
      return false
    end
  end
  return true
end


"""
    $SIGNATURES

Initialization requires down message passing of more specialized down init msgs.
This function performs any possible initialization of variables and retriggers
children cliques that have not yet initialized.

Notes:
- Assumed this function is only called after status from child clique up inits completed.
- Assumes cliq has parent.
  - will fetch message from parent
- Will perform down initialization if status == `:needdownmsg`.
- might be necessary to pass furhter down messges to child cliques that also `:needdownmsg`.
- Will not complete cliq solve unless all children are `:upsolved` (upward is priority).
- `dwinmsgs` assumed to come from parent initialization process.
- assume `subfg` as a subgraph that can be modified by this function (add message factors)
  - should remove message prior factors from subgraph before returning.
- May modify `cliq` values.
  - `setCliqUpInitMsgs!(cliq, cliq.index, msg)`
  - `setCliqStatus!(cliq, status)`
  - `setCliqDrawColor(cliq, "sienna")`
  - `notifyCliqDownInitStatus!(cliq, status)`

Algorithm:
- determine which downward messages influence initialization order
- initialize from singletons to most connected non-singletons
- revert back to needdownmsg if cycleInit does nothing
- can only ever return :initialized or :needdownmsg status
"""
function doCliqInitDown!(subfg::FactorGraph,
                         cliq::Graphs.ExVertex,
                         dwinmsgs::Dict{Symbol,BallTreeDensity})
  #
  @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 1"
  status = :needdownmsg #:badinit
  # get down messages from parent
  # prnt = getParent(tree, cliq)[1]
  # @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 2"
  # dwinmsgs = prepCliqInitMsgsDown!(subfg, tree, prnt)
  @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 3, dwinmsgs=$(collect(keys(dwinmsgs)))"
  # get down variable initialization order
  @show initorder = getCliqInitVarOrderDown(subfg, cliq, dwinmsgs)
  @show map(x->getSym(subfg, x), initorder)

  @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 4, dwinmsgs=$(collect(keys(dwinmsgs)))"
  # add messages as priors to this sub factor graph
  msgfcts = addMsgFactors!(subfg, dwinmsgs)

  @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 5"
  # cycle through vars and attempt init
  if cycleInitByVarOrder!(subfg, initorder)
    # if areCliqVariablesAllInitialized(subfg, cliq)
    status = :initialized
    # end
  end

  @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 6, current status: $status"
  # remove msg factors previously added
  deleteMsgFactors!(subfg, msgfcts)

  @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 7, current status: $status"

  # @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 8, current status: $status"
  # queue the response in a channel to signal waiting tasks
  # @info "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- 9, current status: $status"
  return status
end

function doCliqInitDown!(subfg::FactorGraph,
                         tree::BayesTree,
                         cliq::Graphs.ExVertex  )
  #
  @warn "deprecated doCliqInitDown!(subfg, tree, cliq) use doCliqInitDown!(subfg, cliq, dwinmsgs) instead."
  prnt = getParent(tree, cliq)[1]
  dwinmsgs = prepCliqInitMsgsDown!(subfg, tree, prnt)
  status = doCliqInitDown!(subfg, cliq, dwinmsgs)

  # # TODO move out
  # children = getChildren(tree, cliq)
  # if areCliqChildrenNeedDownMsg(children) # tree, cliq
  #   # status = :initialized
  #   # set messages if children :needdownmsg
  #   @warn "$(current_task()) Clique $(cliq.index), doCliqInitDown! -- must set messages for future down init"
  #   # construct init's up msg to place in parent from initialized separator variables
  #   msg = prepCliqInitMsgsUp(subfg, cliq) # , tree,
  #   @info "$(current_task()) Clique $(cliq.index), putting fake upinitmsg in this cliq, msgs labels $(collect(keys(msg)))"
  #   #fake up message
  #   setCliqUpInitMsgs!(cliq, cliq.index, msg)
  #   setCliqStatus!(cliq, status)
  #   setCliqDrawColor(cliq, "sienna")
  #   notifyCliqDownInitStatus!(cliq, status)
  # end

  return status
end

"""
    $SIGNATURES

Return `true` if any of the children cliques have status `:needdownmsg`.
"""
function areCliqChildrenNeedDownMsg(children::Vector{Graphs.ExVertex})::Bool
  for ch in children
    if getCliqStatus(ch) == :needdownmsg
      return true
    end
  end
  return false
end

function areCliqChildrenNeedDownMsg(tree::BayesTree, cliq::Graphs.ExVertex)::Bool
  areCliqChildrenNeedDownMsg( getChildren(tree, cliq) )
end


"""
    $SIGNATURES

Return true if has parent with status `:needdownmsg`.
"""
function isCliqParentNeedDownMsg(tree::BayesTree, cliq::Graphs.ExVertex)
  prnt = getParent(tree, cliq)
  if length(prnt) == 0
    return false
  end
  prstat = getCliqStatus(prnt[1])
  @info "$(current_task()) Clique $(cliq.index), isCliqParentNeedDownMsg -- parent status: $(prstat)"
  return prstat == :needdownmsg
end

# function doCliqRemoteInitSolve(subfg::FactorGraph)
#
# end

# """
#     $SIGNATURES
#
# Separated function for likely multicore processing, focussed on upward or downward direction initialization of cliques.
#
# Development
# - Make multicore with `remotecall` methods.
# """
# function doCliqInitUpOrDown!(sfg::FactorGraph,
#                              tree::BayesTree,
#                              cliq::Graphs.ExVertex,
#                              isprntnddw::Bool  )
#   #
#   cliqst = getCliqStatus(cliq)
#   # TODO: split if into two states
#   if cliqst == :needdownmsg && !isprntnddw
#     error("obsolete")
#   end
#   if cliqst in [:initialized; :null] && !areCliqChildrenNeedDownMsg(tree, cliq)
#     @info "$(current_task()) Clique $(cliq.index), going for doCliqAutoInitUp!"
#     cliqst = doCliqAutoInitUp!(sfg, tree, cliq)
#   end
#   return (sfg, cliq, cliqst)
# end



# """
#     $SIGNATURES
#
# Major upward initialization / solve inference function.
#
# Notes:
# - will call on values from children or parent cliques
# - can be called multiple times
#
# Future
# - this is a post-hoc (poorly) written state-machine
#   - will rewrite as state machine given benefit of hindsight.
# """
# function cliqInitSolveUp!(fgl::FactorGraph,
#                           tree::BayesTree,
#                           cliq::Graphs.ExVertex;
#                           drawtree::Bool=false,
#                           show::Bool=false,
#                           incremental::Bool=true,
#                           limititers::Int=-1 )
#   #
#   # check clique status
#   cliqst = getCliqStatus(cliq)
#   lbl = cliq.attributes["label"]
#
#   if incremental && cliqst in [:upsolved; :downsolved; :marginalized]
#     # prep and send upward message
#     prnt = getParent(tree, cliq)
#     if length(prnt) > 0
#       # not a root clique
#       # construct init's up msg to place in parent from initialized separator variables
#       msg = prepCliqInitMsgsUp(fgl, tree, cliq)
#       setCliqUpInitMsgs!(prnt[1], cliq.index, msg)
#       notifyCliqUpInitStatus!(cliq, cliqst)
#       @info "$(current_task()) Clique $(cliq.index), skip computation on status=$cliqst, but did prepare/notify upward message"
#     end
#
#     return cliqst
#   end
#
#   # build a local subgraph for inference operations
#   syms = getCliqAllVarSyms(fgl, cliq)
#   sfg = buildSubgraphFromLabels(fgl, syms)
#
#   # get parent cliq
#   prnt = getParent(tree, cliq)
#
#   tryonce = true
#   countiters = 0
#   # upsolve delay loop
#   while (0 < limititers || limititers == -1) && (tryonce || !(cliqst in [:upsolved; :downsolved; :marginalized]))
#     countiters += 1
#     @info "$(current_task()) Clique $(cliq.index), #$countiters, top of while"
#     limititers != -1 ? (limititers -= 1) : nothing
#     tryonce = false
#     forceproceed = false
#     cliqst = getCliqStatus(cliq)
#     stdictprnt = Dict{Int, Symbol}()
#
#     # @info "$(current_task()) Clique $(cliq.index), status $cliqst -- top of while loop"
#     if cliqst == :needdownmsg && length(prnt) > 0
#       # wait here until all children have a valid status
#       if !areCliqChildrenNeedDownMsg(tree, cliq)
#         @info "$(current_task()) Clique $(cliq.index), blocking on parent until all sibling cliques have valid status"
#         setCliqDrawColor(cliq, "turquoise")
#         drawtree ? drawTree(tree, show=show) : nothing
#         stdictprnt = blockCliqUntilChildrenHaveUpStatus(tree, prnt[1])
#       else
#         @warn "$(current_task()) Clique $(cliq.index), WIP must deal with child :needdownmsg"
#         forceproceed = true
#       end
#     end
#
#     # Determine if child clique processes all completed with status :upsolved
#     @info "$(current_task()) Clique $(cliq.index), cliqInitSolveUp! -- blocking until child cliques have status, cliqst=$(cliqst)"
#     stdict = blockCliqUntilChildrenHaveUpStatus(tree, cliq)
#     @info "$(current_task()) Clique $(cliq.index) continue, children all have status"
#
#     # promote if longer down chain of :needdownmsg
#     if cliqst == :null
#       chstatus = collect(values(stdict))
#       len = length(chstatus)
#       if len > 0 && sum(chstatus .== :needdownmsg) == len
#         # TODO maybe can happen where some children need more information?
#         @info "$(current_task()) Clique $(cliq.index) | $lbl | escalating to :needdownmsg since all children :needdownmsg"
#         notifyCliqUpInitStatus!(cliq, :needdownmsg)
#         # setCliqStatus!(cliq, :needdownmsg)
#         cliqst = getCliqStatus(cliq)
#         setCliqDrawColor(cliq, "green")
#         tryonce = true
#       end
#
#       # wait if child branches still solving -- must eventually upsolve this clique
#       if len > 0 && sum(chstatus .!= :upsolved) > 0
#         @info "$(current_task()) Clique $(cliq.index) | $lbl | sleeping until all children finish upward inference"
#         sleep(0.1)
#       end
#     end
#
#     # hard assumption here on upsolve from leaves to root
#     proceed = true
#     # TODO not sure if we want stdict from cliq or prnt???
#     for (clid, clst) in stdict
#       @info "$(current_task()) Clique $(cliq.index), check stdict: clid=$(clid), clst=$(clst)"
#       # :needdownmsg # 'send' downward init msg direction
#       # :initialized # @warn "something might not be right with init of clid=$clid"
#       !(clst in [:initialized;:upsolved;:marginalized;:downsolved]) ? (proceed = false) : nothing
#     end
#     @info "$(current_task()) Clique $(cliq.index), proceed=$(proceed), tryonce=$tryonce, clst=$(cliqst)"
#
#     # add blocking case when all siblings and parent :needdownmsg -- until parent :initialized
#     @info "$(current_task()) Clique $(cliq.index), check block if siblings & parent have :needdownmsg status? clst=$(cliqst), proceed=$proceed, forceproceed=$forceproceed."
#     blockCliqSiblingsParentNeedDown(tree, cliq)
#
#     # add case for if children are blocked on need down msg
#     if getCliqStatus(cliq) == :initialized && areCliqChildrenNeedDownMsg(tree, cliq)
#       sleep(0.1)
#     end
#
#     # if all children are ready, proceed with this cliq initialization
#     if proceed || forceproceed
#       # start computations
#       setCliqDrawColor(cliq, "red")
#       drawtree ? drawTree(tree, show=show) : nothing
#       # evaluate according to cliq status
#       isprntnddw = isCliqParentNeedDownMsg(tree, cliq)
#       @info "$(current_task()) Clique $(cliq.index), proceed: $(cliqst), isCliqParentNeedDownMsg(tree, cliq)=$(isprntnddw), areCliqChildrenNeedDownMsg(tree, cliq)=$(areCliqChildrenNeedDownMsg(tree, cliq))"
#       d1,d2,cliqst = doCliqInitUpOrDown!(sfg, tree, cliq, isprntnddw)
#       # if cliqst == :needdownmsg && !isprntnddw
#       #   # initialize clique in downward direction
#       #   # not if parent also needs downward init message
#       #   @info "$(current_task()) Clique $(cliq.index), needs down message -- attempt down init"
#       #   cliqst = doCliqInitDown!(sfg, tree, cliq)
#       #   @info "$(current_task()) Clique $(cliq.index), after down init attempt, $cliqst."
#       # end
#       # if cliqst in [:initialized; :null] && !areCliqChildrenNeedDownMsg(tree, cliq)
#       #   @info "$(current_task()) Clique $(cliq.index), going for doCliqAutoInitUp!"
#       #   cliqst = doCliqAutoInitUp!(sfg, tree, cliq)
#       # end
#       if cliqst == :upsolved
#         @info "$(current_task()) Clique $(cliq.index), going for transferUpdateSubGraph!"
#         frsyms = Symbol[getSym(sfg, varid) for varid in getCliqFrontalVarIds(cliq)]
#         transferUpdateSubGraph!(fgl, sfg, frsyms)
#       elseif cliqst == :initialized
#         # @info "$(current_task()) Clique $(cliq.index), set update down init messages: "  # OBSOLETE
#         setCliqDrawColor(cliq, "sienna")
#       else
#         @info "$(current_task()) Clique $(cliq.index), init not complete and should wait on init down message."
#         setCliqDrawColor(cliq, "green")
#         tryonce = true # TODO, potential problem with trying to downsolve
#       end
#       drawtree ? drawTree(tree, show=show) : nothing
#     end
#     @info "$(current_task()) Clique $(cliq.index), #$countiters, bottom of while, cliqst=$(cliqst)"
#   end # while
#   @info "$(current_task()) Clique $(cliq.index), total #$countiters, after while completed up inference."
#   return cliqst
# end
#
#
# """
#     $SIGNATURES
#
# Based on current status in factor graph, determine if initialization of requested
# variable is possible.
# """
# function getCliqInitVarPossible(cliq::Graphs.ExVertex, varid::Int)
#
#   factorCanInitFromOtherVars(cliq, fctid)
#
# end

# """
#     $SIGNATURES
#
# Determine if a clique Chapman-Kolmogorov computation can be achieved,
# alongide additional message singletons that might be available from caller.
# """
# function calcCliqTotalSolvePossible(cliq::Graphs.ExVertex;
#                                     allids::Vector{Int}=getCliqAllVarIds(cliq),
#                                     availablemsgs::Vector{Bool}=zeros(Bool,length(allids)) )::Tuple{Bool, Vector{Bool}}
#   # return list of all initable variables in cliq (default is false)
#   initable = zeros(Bool, length(allids))
#
#   # what is the initialization order
#   initorder = getCliqInitVarOrderUp(cliq::Graphs.ExVertex)
#
#   # check if all variables can be initialized
#   for i in 1:length(allids)
#     if allids[i] in initorder
#       if getCliqInitVarPossible(cliq, allids[i])
#         initable[i] = true
#       end
#     end
#   end
#
#   # would be fully initializable if all initable are true
#   return sum(initable)==length(allids), initable
# end
