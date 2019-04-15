
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

getCliqInitDownMsgs(cliq::Graphs.ExVertex) = getData(cliq).downInitMsg


"""
    $SIGNATURES

Return the most likely  ordering for initializing factor (assuming up solve
sequence).
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

Perform cliq initalization calculation based on current state of the tree and factor graph,
using upward message passing logic.

> NOTE WORK IN PROGRESS

Notes
- Return either of (:initialized, :upsolved, :needdownmsg, :badinit)
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
  count = 1
  while count > 0
    count = 0
    for vid in varorder
      var = getVert(fgl, vid, api=localapi)
      isinit = isInitialized(var)
      doautoinit!(fgl, ExVertex[var;], api=localapi)
      isinit == isInitialized(var) ? nothing : (count += 1)
    end
  end

  # next default return type
  retmsg = :needdownmsg

  # check if all cliq vars have been initialized so that full inference can occur on clique
  isinit = areCliqVariablesInitialized(fgl, cliq)
  # might fail while waiting for other cliques to initialize.
  if isinit
    retmsg = :initialized
    if up_solve_if_able
      csym = Symbol(getVert(fgl,getCliqFrontalVarIds(cliq)[1],api=localapi).label)
      approxCliqMarginalUp!(fgl, tree, csym, false)
      retmsg = :upsolved
    end
  end

  # construct init msg to place in parent as initialized separator variables
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

Block the thread until child cliques of `prnt::Graphs.ExVertex` have been initialized.
Return `::Symbol` indicating whether next action that should be taken, namely:
- :needdownmsg
- :ready_upsolve
- :ready_initialization
- :unknown

Notes:
- Can be called multiple times
"""
function blockUntilCliqChildrenUpInit(tree::BayesTree, prnt::Graphs.ExVertex)
  retmsg = :unknown
  chlr = getChildren(tree, prnt)
  for ch in chlr
    if !isready(getData(ch).initUpChannel) && isCliqInitialized(ch)

    end

  end

  return
end

"""
    $SIGNATURES

Initialization downward message passing is different from regular inference since
it is possible that none of the child cliq variables have been initialized.
"""
function prepCliqInitMsgsDown!(fgl::FactorGraph, tree::BayesTree, cliq::Graphs.ExVertex)
  #
  # get the parent cliq
  prnt = getParent(tree, cliq)
  # get the current messages stored in the parent
  currmsgs = getCliqInitUpMsgs(pcliq)

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

function doCliqAutoInitDown!(fgl::FactorGraph,
                             tree::BayesTree,
                             cliq::Graphs.ExVertex   )
  #
  status = :starting
  # get down messages from parent
  prnt = getParent(tree, cliq)
  msgs = getCliqInitDownMsgs(prnt)
  cliqd = getData(cliq)


  @warn "work in progress"

  put!(cliqd.initDownChannel, status)
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
