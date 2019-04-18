
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

"""
    $SIGNATURES

Return true if clique has completed the local upward direction inference procedure.
"""
isUpInferenceComplete(cliq::Graphs.ExVertex) = getData(cliq).upsolved

function areCliqVariablesInitialized(fgl::FactorGraph, cliq::Graphs.ExVertex)
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
  isallinit = areCliqVariablesInitialized(fgl, cliq)

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

Cycle through var order and initialize variables as possible in `subfg::FactorGraph`.

Notes:
- assumed `subfg` is a subgraph containing only the factors that can be used.
- intended for both up and down initialization operations.
"""
function cycleInitByVarOrder!(subfg::FactorGraph, varorder::Vector{Int})
  count = 1
  while count > 0
    count = 0
    for vid in varorder
      var = getVert(subfg, vid, api=localapi)
      isinit = isInitialized(var)
      # TODO -- must use factors and values in cliq (assume subgraph?)
      doautoinit!(subfg, ExVertex[var;], api=localapi)
      isinit == isInitialized(var) ? nothing : (count += 1)
    end
  end
  nothing
end

"""
    $SIGNATURES

Return (:badinit, :initialized, upsolved) whether a cliq has been upsolved
(possibly even before this call), and update `subfg::FactorGraph` according to
internal computations from this function.
"""
function attemptCliqUpSolve!(subfg::FactorGraph,
                             tree::BayesTree,
                             cliq::Graphs.ExVertex)::Symbol
  #
  retmsg = :needdownmsg
  # check if all cliq vars have been initialized so that full inference can occur on clique
  isinit = areCliqVariablesInitialized(subfg, cliq)
  # might fail while waiting for other cliques to initialize.
  if isinit
    retmsg = :initialized
    csym = Symbol(getVert(subfg, getCliqFrontalVarIds(cliq)[1], api=localapi).label)
    approxCliqMarginalUp!(subfg, tree, csym, false)
    retmsg = :upsolved
  end
  return retmsg
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
function doCliqAutoInitUp!(fgl::FactorGraph,
                           tree::BayesTree,
                           cliq::Graphs.ExVertex;
                           up_solve_if_able::Bool=true  )::Symbol
  #
  # init up msg has special procedure for incomplete messages
  retmsg = :badinit
  msg = Dict{Symbol, BallTreeDensity}()

  # structure for all up message densities computed during this initialization procedure.
  varorder = getCliqInitVarOrderUp(cliq)

  # do physical inits
  cycleInitByVarOrder!(fgl, varorder)

  # # check if all cliq vars have been initialized so that full inference can occur on clique
  retmsg = attemptCliqUpSolve!(fgl, tree, cliq)

  # construct init's up msg to place in parent from initialized separator variables
  for vid in getCliqSeparatorVarIds(cliq)
    var = getVert(fgl, vid, api=localapi)
    if isInitialized(var)
      msg[Symbol(var.label)] = getKDE(var)
    end
  end

  # put the init result in the parent cliq.
  prnt = getParent(tree, cliq)
  # not a root clique
  if length(prnt) > 0
    setCliqUpInitMsgs!(prnt[1], cliq.index, msg)
  end

  # set flags in clique for multicore sequencing
  getData(cliq).initialized = retmsg
  put!(getData(cliq).initUpChannel, retmsg)
  return retmsg
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
  ret = Dict{Int, Symbol}()
  chlr = getChildren(tree, prnt)
  for ch in chlr
    # either wait to fetch new result, or report or result
    if !isready(getData(ch).initUpChannel) && getCliqStatusUp(ch) != :null
      ret[ch.index] = getCliqStatusUp(ch)
    else
      ret[ch.index] = take!(getData(ch).initUpChannel)
    end
  end
  return ret
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
    if length(msgspervar[msgsym]) > 1
      products[msgsym] = manifoldProduct(msgs, getManifolds(fgl, msgsym))
    else
      products[msgsym] = msgs[1]
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

Initialization requires down message passing of more specialized down init msgs.
This function performs any possible initialization of variables and retriggers
children cliques that have not yet initialized.

Notes:
- Assumed this function is only called after status from child clique up inits completed.
- Will perform down initialization if status == `:needdownmsg`.
- might be necessary to pass furhter down messges to child cliques that also `:needdownmsg`.
- Will not complete cliq solve unless all children are `:upsolved` (upward is priority).
- `dwinmsgs` assumed to come from parent initialization process.

Algorithm:
- determine which downward messages influence initialization order
"""
function doCliqAutoInitDown!(fgl::FactorGraph,
                             tree::BayesTree,
                             cliq::Graphs.ExVertex  )
  #
  status = :badinit
  # get down messages from parent
  @show prnt = getParent(tree, cliq)[1]
  dwinmsgs = prepCliqInitMsgsDown!(fgl, tree, prnt)

  # get down variable initialization order
  @show initorder = getCliqInitVarOrderDown(fgl, cliq, dwinmsgs)

  # cycle through vars and attempt init
  cycleInitByVarOrder!(fgl, initorder)

  # check if cliq variables have been initialized
  @show tempstatus = areCliqVariablesInitialized(fgl, cliq)

  # check if subset of children have initialized or solved.


  # check if any child cliques `:needdownmsg`, and prepare accordingly

  @warn "work in progress"

  # put!(cliqd.initDownChannel, status)
  return status
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
