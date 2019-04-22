
"""
    $SIGNATURES

Based on a push model from child cliques that should have already completed their computation.
"""
function getCliqInitUpMsgs(cliq::Graphs.ExVertex)
  getData(cliq).upInitMsgs
end

function setCliqUpInitMsgs!(cliq::Graphs.ExVertex, childid::Int, msg::Dict{})
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
    @info "dumping stale cliq=$(cliq.index) status message $(take!(cd.initUpChannel)), replacing with $(status)"
  end
  put!(cd.initUpChannel, status)
end

function notifyCliqDownInitStatus!(cliq::Graphs.ExVertex, status::Symbol)
  cd = getData(cliq)
  cd.initialized = status
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
    @show isallinit &= isUpInferenceComplete(chl)
  end
  return isallinit
end

"""
    $SIGNATURES

Blocking call until `cliq` upInit processes has arrived at a result.
"""
function getCliqInitUpResultFromChannel(cliq::Graphs.ExVertex)
  take!(getData(cliq).initUpChannel)
end

"""
    $SIGNATURES

Return `::Symbol` status a particular clique is in, with specific regard to solution
or numerical initialization status:
- :needdownmsg
- :upsolved
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

Block the thread until child cliques of `prnt::Graphs.ExVertex` have finished
attempting upward initialization -- i.e. have status result.
Return `::Dict{Symbol}` indicating whether next action that should be taken
for each child clique.

Notes:
- See status options at `getCliqStatusUp(..)`.
- Can be called multiple times
"""
function blockCliqUntilChildrenHaveUpStatus(tree::BayesTree,
                                            prnt::Graphs.ExVertex)::Dict{Int, Symbol}
  #
  @warn "in blockCliqUntilChildrenHaveUpStatus"
  ret = Dict{Int, Symbol}()
  chlr = getChildren(tree, prnt)
  for ch in chlr
    # either wait to fetch new result, or report or result
    if !isready(getData(ch).initUpChannel) && getCliqStatusUp(ch) != :null
      ret[ch.index] = getCliqStatusUp(ch)
    else
      ret[ch.index] = fetch(getData(ch).initUpChannel) # take!
    end
  end
  return ret
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
  retval = false
  count = 1
  while count > 0
    count = 0
    for vid in varorder
      var = getVert(subfg, vid, api=localapi)
      isinit = isInitialized(var)
      # TODO -- must use factors and values in cliq (assume subgraph?)
      doautoinit!(subfg, ExVertex[var;], api=localapi)
      if isinit != isInitialized(var)
        count += 1
        retval = true
      end
    end
  end
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
function prepCliqInitMsgsUp!(subfg::FactorGraph, tree::BayesTree, cliq::Graphs.ExVertex)::Dict{Symbol, BallTreeDensity}
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
  @info "adding up message factors"
  msgfcts = addMsgFactors!(subfg, upmsgs)

    # TEMP
    writeGraphPdf(subfg, show=true)

  # attempt initialize if necessary
  if !areCliqVariablesAllInitialized(subfg, cliq)
    # structure for all up message densities computed during this initialization procedure.
    varorder = getCliqInitVarOrderUp(cliq)
    # do physical inits, ignore cycle return value
    @show "going for up cycle order"
    cycleInitByVarOrder!(subfg, varorder)
    @show "finished with up cycle order"
  end

  # check if all cliq vars have been initialized so that full inference can occur on clique
  if areCliqVariablesAllInitialized(subfg, cliq)
    @show status = doCliqUpSolve!(subfg, tree, cliq)
  end

  # construct init's up msg to place in parent from initialized separator variables
  msg = prepCliqInitMsgsUp!(subfg, tree, cliq)

  # put the init result in the parent cliq.
  @show prnt = getParent(tree, cliq)
  if length(prnt) > 0
    # not a root clique
    setCliqUpInitMsgs!(prnt[1], cliq.index, msg)
  end

  # remove msg factors that were added to the subfg
  @info "removing up message factors, length=$(length(msgfcts))"
  deleteMsgFactors!(subfg, msgfcts)

  # set flags in clique for multicore sequencing
  @info "sending notification of up init status"
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
function prepCliqInitMsgsDown!(fgl::FactorGraph, tree::BayesTree, cliq::Graphs.ExVertex)
  #
  # get the current messages stored in the parent
  currmsgs = getCliqInitUpMsgs(cliq)

  # check if any msgs should be multiplied together for the same variable
  msgspervar = Dict{Symbol, Vector{BallTreeDensity}}()
  for (cliqid, msgs) in currmsgs
    for (msgsym, msg) in msgs
      if !haskey(msgspervar, msgsym)
        msgspervar[msgsym] = Vector{BallTreeDensity}()
      end
      push!(msgspervar[msgsym], msg)
    end
  end

  # multiply multiple messages together
  products = getData(cliq).downInitMsg # Dict{Symbol, BallTreeDensity}()
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
- Will perform down initialization if status == `:needdownmsg`.
  - will fetch message from parent
- might be necessary to pass furhter down messges to child cliques that also `:needdownmsg`.
- Will not complete cliq solve unless all children are `:upsolved` (upward is priority).
- `dwinmsgs` assumed to come from parent initialization process.
- assume `subfg` as a subgraph that can be modified by this function (add message factors)
  - should remove message factors from subgraph before returning.

Algorithm:
- determine which downward messages influence initialization order
- initialize from singletons to most connected non-singletons
"""
function doCliqInitDown!(subfg::FactorGraph,
                         tree::BayesTree,
                         cliq::Graphs.ExVertex  )
  #
  @warn "entering doCliqInitDown"
  status = :badinit
  # get down messages from parent
  @show prnt = getParent(tree, cliq)[1]
  dwinmsgs = prepCliqInitMsgsDown!(subfg, tree, prnt)

  # get down variable initialization order
  @show initorder = getCliqInitVarOrderDown(subfg, cliq, dwinmsgs)

  # add messages as priors to this sub factor graph
  msgfcts = addMsgFactors!(subfg, dwinmsgs)

  # cycle through vars and attempt init
  if cycleInitByVarOrder!(subfg, initorder)
    @show status = :initialized
  end

  # check if all cliq variables have been initialized
  if !areCliqVariablesAllInitialized(subfg, cliq)
    @show status = :needdownmsg
  end

  # remove msg factors previously added
  deleteMsgFactors!(subfg, msgfcts)

  # queue the response in a channel to signal waiting tasks
  notifyCliqDownInitStatus!(cliq, status)
  return status
end


"""
    $SIGNATURES

Major upward initialization / solve inference function.

Notes:
- will call on values from children or parent cliques
- can be called multiple times
"""
function cliqInitSolveUp!(fgl::FactorGraph,
                          tree::BayesTree,
                          cliq::Graphs.ExVertex;
                          drawtree::Bool=false,
                          show::Bool=false,
                          incremental::Bool=true  )
  #
  # check clique status
  cliqst = getCliqStatus(cliq)

  if incremenal && cliqst in [:upsolved; :downsolved; :marginalized]
    return cliqst
  end

  syms = getCliqAllVarSyms(fgl, cliq)
  # build a local subgraph for inference operations
  sfg = buildSubgraphFromLabels(fgl, syms)

  # get parent cliq
  prnt = getParent(tree, cliq)

  tryonce = true
  # upsolve delay loop
  while tryonce || !areCliqChildrenAllUpSolved(tree, cliq)
    tryonce = false
    cliqst = getCliqStatus(cliq)

    # @info "clique $(cliq.index), status $cliqst -- top of while loop"
    if cliqst == :needdownmsg && length(prnt) > 0
      # wait here until all children have a valid status
      @info "clique $(cliq.index), blocking on parent until sibling cliques have valid status"
      setCliqDrawColor(cliq, "turquoise")
      drawtree ? drawTree(tree, show=show) : nothing
      stdict = blockCliqUntilChildrenHaveUpStatus(tree, prnt[1])
      # @info "clique $(cliq.index), siblings have status"
    end

    # Determine if child clique processes all completed with status :upsolved
    @info "blocking on clique $(cliq.index) until child cliques have status"
    stdict = blockCliqUntilChildrenHaveUpStatus(tree, cliq)
    # @info "clique $(cliq.index) continue, children all have status"

    # hard assumption here on upsolve from leaves to root
    proceed = true
    for (clid, clst) in stdict
      # :needdownmsg # 'send' downward init msg direction
      # :initialized # @warn "something might not be right with init of clid=$clid"
      !(clst in [:initialized;:upsolved]) ? (proceed = false) : nothing
    end

    # if all children are ready, proceed with this cliq initialization
    if proceed
      # start computations
      setCliqDrawColor(cliq, "red")
      drawtree ? drawTree(tree, show=show) : nothing
      # evaluate according to cliq status
      if cliqst == :needdownmsg
        # initialize clique in downward direction
        @info "this clique needs down message information and will now try do down init."
        cliqst = doCliqInitDown!(sfg, tree, cliq)
      end
      if cliqst in [:initialized; :null]
        @info "going for doCliqAutoInitUp!"
        cliqst = doCliqAutoInitUp!(sfg, tree, cliq)
      end
      if cliqst == :upsolved
        @info "going for transferUpdateSubGraph!"
        frsyms = Symbol[getSym(sfg, varid) for varid in getCliqFrontalVarIds(cliq)]
        transferUpdateSubGraph!(fgl, sfg, frsyms)
      else
        @info "clique $(cliq.index), init waiting since it cannot fully up solve yet."
        setCliqDrawColor(cliq, "green")
        tryonce = true
      end
      drawtree ? drawTree(tree, show=show) : nothing
    end
  end # while
  return cliqst
end



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
