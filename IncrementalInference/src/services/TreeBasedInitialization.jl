
function isCliqInitialized(cliq::TreeClique)::Bool
  return getCliqueData(cliq).initialized in [INITIALIZED; UPSOLVED]
end

function isCliqUpSolved(cliq::TreeClique)::Bool
  return getCliqueData(cliq).initialized == UPSOLVED
end

"""
    $SIGNATURES

Return the most likely  ordering for initializing factor (assuming up solve
sequence).

Notes:
- sorts id (label) for increasing number of connected factors using the clique subfg with messages already included.
"""
function getCliqVarInitOrderUp(subfg::AbstractDFG)
  # rules to explore dimension from one to the other?

  # get all variable ids and number of associated factors
  B, varLabels, facLabels = getBiadjacencyMatrix(subfg)
  nfcts = sum(B; dims = 1)[:]

  # variables with priors
  varswithpriors = listNeighbors.(subfg, lsfPriors(subfg))
  singids = union(Symbol[], varswithpriors...)

  # sort permutation order for increasing number of factor association
  nfctsp = sortperm(nfcts)
  sortedids = varLabels[nfctsp]

  # organize the prior variables separately with asceding factor count
  initorder = intersect(sortedids, singids)
  # in ascending order of number of factors
  union!(initorder, sortedids)

  return initorder
end

"""
    $SIGNATURES

Special function to do initialization in downward direction, assuming that not all
variables can be initialized.  Relies on outside down messages.

Notes:
- assumed this `cliq` is being initialized from a previous `:needdownmsg` status.
- will use all possible local factors of cliquq in initilization process
- similar to upward initialization, but uses different message structure
  - first draft assumes upward messages will not be used,
  - full up solve still required which explicitly depends on upward messages.

Dev Notes
- Streamline add/delete msg priors from calling function and csm.
- TODO replace with nested 'minimum degree' type variable ordering.
"""
function getCliqInitVarOrderDown(
  dfg::AbstractDFG,
  cliq::TreeClique,
  dwnkeys::Vector{Symbol},
)   # downmsgs
  #
  allsyms = getCliqAllVarIds(cliq)
  # convert input downmsg var symbols to integers (also assumed as prior beliefs)
  # make sure ids are in the clique set, since parent may have more variables.
  dwnmsgsym = intersect(dwnkeys, DFG.listVariables(dfg))
  dwnvarids = intersect(allsyms, dwnmsgsym)

  # find any other prior factors (might have partials)
  prvarids = getCliqVarIdsPriors(cliq, allsyms, true)
  hassinglids = union(dwnvarids, prvarids)

  # Get all other variable factor counts
  nfcts = getCliqNumAssocFactorsPerVar(cliq)
  # add msg marginal prior (singletons) to number of factors
  for msid in dwnmsgsym
    nfcts[msid .== allsyms] .+= 1
  end

  # sort permutation order for increasing number of factor association
  nfctsp = sortperm(nfcts)
  sortedids = allsyms[nfctsp]

  # all singleton variables
  singids = union(prvarids, dwnvarids)

  # organize the prior variables separately with asceding factor count
  initorder = Symbol[] #zeros(Int, 0)
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
  return initorder::Vector{Symbol}
end

function _isInitializedOrInitSolveKey(
  var::VariableCompute,
  solveKey::Symbol = :default;
  N::Int = 100,
)
  # TODO, this solveKey existence test should probably be removed?
  if !(solveKey in listStates(var))
    varType = getStateKind(var)
    setDefaultNodeData!(
      var,
      0,
      N;
      solveKey = solveKey,
      initialized = false,
      varType = varType,
      # dontmargin = false,
    )
    #
    # data = getState(var, solveKey)
    # if data === nothing
    # end
    return false
  end
  # regular is initialized check, this is fine
  isinit = isInitialized(var, solveKey)
  return isinit
end

"""
    $SIGNATURES

Return true if clique has completed the local upward direction inference procedure.
"""
isUpInferenceComplete(cliq::TreeClique) = getCliqueData(cliq).upsolved

function areCliqVariablesAllInitialized(
  dfg::AbstractDFG,
  cliq::TreeClique,
  solveKey::Symbol = :default;
  N::Int = getSolverParams(dfg).N,
)
  #
  allids = getCliqAllVarIds(cliq)
  isallinit = true
  for vid in allids
    var = DFG.getVariable(dfg, vid)
    isallinit &= _isInitializedOrInitSolveKey(var, solveKey; N = N)
    # isallinit &= isInitialized(var, solveKey)
  end
  return isallinit
end

"""
    $SIGNATURES

Return true if all variables in clique are considered marginalized (and initialized).
"""
function areCliqVariablesAllMarginalized(subfg::AbstractDFG, cliq::TreeClique)
  for vsym in getCliqAllVarIds(cliq)
    vert = getVariable(subfg, vsym)
    if !isMarginalized(vert) || !isInitialized(vert)
      return false
    end
  end
  return true
end

function printCliqInitPartialInfo(
  subfg,
  cliq,
  solveKey::Symbol = :default,
  logger = ConsoleLogger(),
)
  varids = getCliqAllVarIds(cliq)
  initstatus = Vector{Bool}(undef, length(varids))
  initpartial = Vector{Float64}(undef, length(varids))
  for i = 1:length(varids)
    initstatus[i] = isInitialized(subfg, varids[i], solveKey) # getState(getVariable(subfg, varids[i]), solveKey).initialized
    initpartial[i] = -1 # getState(getVariable(subfg, varids[i]), solveKey).inferdim
  end
  with_logger(logger) do
    tt = split(string(now()), 'T')[end]
    @info "$tt, cliq $(cliq.id), PARINIT: $varids | $initstatus | $initpartial"
  end
end

#
