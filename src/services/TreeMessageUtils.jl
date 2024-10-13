# init utils for tree based inference

## =============================================================================
# short preamble funcions
## =============================================================================

function convert(::Type{<:ManifoldKernelDensity}, src::TreeBelief)
  return manikde!(getManifold(src.variableType), src.val; bw = src.bw[:, 1])
end

manikde!(em::TreeBelief) = convert(ManifoldKernelDensity, em)

## =============================================================================
# helper functions for tree message channels
## =============================================================================

"""
    $SIGNATURES

Reset the state of all variables in a clique to not initialized.

Notes
- resets numberical values to zeros.

Dev Notes
- TODO not all kde manifolds will initialize to zero.
- FIXME channels need to be consolidated
"""
function resetCliqSolve!(
  dfg::AbstractDFG,
  treel::AbstractBayesTree,
  cliq::TreeClique;
  solveKey::Symbol = :default,
)
  #
  cda = getCliqueData(cliq)
  vars = getCliqVarIdsAll(cliq)
  for varis in vars
    resetVariable!(dfg, varis; solveKey = solveKey)
  end
  # TODO remove once consolidation with upMsgs is done
  putCliqueMsgUp!(cda, LikelihoodMessage())

  # cda.dwnMsg = LikelihoodMessage()
  putCliqueInitMsgDown!(cda, LikelihoodMessage())

  setCliqueStatus!(cliq, NULL)
  setCliqueDrawColor!(cliq, "")
  return nothing
end

function resetCliqSolve!(
  dfg::AbstractDFG,
  treel::AbstractBayesTree,
  frt::Symbol;
  solveKey::Symbol = :default,
)
  #
  return resetCliqSolve!(dfg, treel, getClique(treel, frt); solveKey = solveKey)
end

## =============================================================================
# helper functions to add tree messages to subgraphs
## =============================================================================

function updateSubFgFromDownMsgs!(
  sfg::G,
  dwnmsgs::LikelihoodMessage,
  seps::Vector{Symbol},
) where {G <: AbstractDFG}
  #
  # sanity check basic Bayes (Junction) tree property
  # length(setdiff(keys(dwnmsgs), seps)) == 0 ? nothing : error("updateSubFgFromDownMsgs! -- separators and dwnmsgs not consistent")

  # update specific variables in sfg from msgs
  for (key, beldim) in dwnmsgs.belief
    if key in seps
      newBel = manikde!(getManifold(beldim.variableType), beldim.val; bw = beldim.bw[:, 1])
      setValKDE!(sfg, key, newBel, false, beldim.infoPerCoord)
    end
  end

  return nothing
end

function generateMsgPrior(belief_::TreeBelief, ::NonparametricMessage)
  kdePr = manikde!(getManifold(belief_.variableType), belief_.val; bw = belief_.bw[:, 1])
  return MsgPrior(kdePr, belief_.infoPerCoord, getManifold(belief_))
end

function generateMsgPrior(belief_::TreeBelief, ::ParametricMessage)
  msgPrior = if length(belief_.val[1]) == 1 #FIXME ? && length(belief_.val) == 1
    MsgPrior(
      Normal(belief_.val[1][1], sqrt(belief_.bw[1])),
      belief_.infoPerCoord,
      getManifold(belief_),
    )
  elseif length(belief_.val[1]) > 1 #FIXME ? length(belief_.val) == 1
    mvnorm = createMvNormal(belief_.val[1], belief_.bw)
    mvnorm !== nothing ? nothing : (return DFGFactor[])
    MsgPrior(mvnorm, belief_.infoPerCoord, getManifold(belief_))
  end
  return msgPrior
end

"""
    $SIGNATURES

Return `Dict{Int, Vector{Symbol}}` where each `Int` is a new subgraph and the vector contains all variables
connected to that subgraph.  Subgraphs connectivity is defined by factors of the [`selectFactorType`](@ref) 
type -- e.g. `Pose2` variables connected by a chain of `Pose2Pose2` factors is connected, but not if a link
is `Pose2Pose2Range`.  This function is specifically intended for use with `MessageRelativeLikelihoods` in mind
to determine which relative and prior factors should be included in an upward belief propagation (joint) message.
Each returned subgraph should receive a `MsgPrior` on the dominant variable.

Notes
- Disconnected subgraphs in the separator variables of a clique should instead be connected by a 
  `TangentAtlasFactor` approximation -- i.e. when analytical `selectFactorType`s cannot be used.
- Internally does `getfield(Main, Symbol(factorname::Core.TypeName))` which might cause unique situations with weird user functions
  - As well as a possible speed penalty -- TODO, investigate

Related

[`_calcCandidatePriorBest`](@ref)
"""
function _findSubgraphsFactorType(
  dfg_::AbstractDFG,
  jointrelatives::MsgRelativeType,
  separators::Vector{Symbol},
)
  #
  commonJoints = []
  subClassify = Dict{Symbol, Int}()
  newClass = 0

  # 1. count separtor connectivity in UPWARD_DIFFERENTIAL
  sepsCount = Dict{Symbol, Int}()
  map(x -> (sepsCount[x] = 0), separators)
  # tagsFilter = [:__LIKELIHOODMESSAGE__;]
  # tflsf = lsf(fg, tags=tagsFilter)
  for likl in jointrelatives
    for vari in likl.variables
      sepsCount[vari] += 1
    end
  end

  # 2. start with 0's as subgraphs
  for (id, count) in sepsCount
    if count == 0
      # also keep second list just to be sure based on labels
      newClass += 1
      subClassify[id] = newClass
    end
  end

  # 3. then < 0 and search all paths, adding each hit to subgraph classifications
  for key1 in setdiff(keys(sepsCount), keys(subClassify))
    if !(key1 in keys(subClassify))
      newClass += 1
      subClassify[key1] = newClass
    end
    # if sepsCount[key1] == 1
    # search connectivity throughout remaining variables, some duplicate computation occurring
    for key2 in setdiff(keys(sepsCount), keys(subClassify))
      defaultFct = selectFactorType(dfg_, key1, key2)
      # @show key1, key2, defaultFct
      # TODO validate getfield Main here
      # resname = defaultFct isa UnionAll ? getfield(Main, defaultFct.body.name |> Symbol) : defaultFct
      resname =
        defaultFct isa UnionAll ? getfield(Main, nameof(defaultFct.body)) : defaultFct
      pth = findShortestPathDijkstra(
        dfg_,
        key1,
        key2;
        typeFactors = [resname;],
        initialized = true,
      )
      # check if connected to existing subClass
      if 0 == length(pth)
        # not connected, so need new class
        newClass += 1
        subClassify[key2] = newClass
      else
        # is connected, so add existing class of key1
        subClassify[key2] = subClassify[key1]
      end
    end
    # end
  end

  # 4. inverse classification dictionary
  allClasses = Dict{Int, Vector{Symbol}}()
  for (key, cls) in subClassify
    if isInitialized(dfg_, key)
      if !haskey(allClasses, cls)
        (allClasses[cls] = Symbol[key;])
      else
        union!(allClasses[cls], [key;])
      end
    end
  end

  # 
  return allClasses
end

"""
    $SIGNATURES

Build from a `LikelihoodMessage` a temporary distributed factor graph object containing differential
information likelihood factors based on values in the messages.

Notes
- Modifies tfg argument by adding `:__UPWARD_DIFFERENTIAL__` factors.

DevNotes
- Initial version which only works for Pose2 and Point2 at this stage.
"""
function addLikelihoodsDifferential!(
  msgs::LikelihoodMessage,
  cliqSubFG::AbstractDFG,
  tfg::AbstractDFG = initfg(),
)
  # create new local dfg and add all the variables with data

  for difflikl in msgs.jointmsg.relatives
    addFactor!(
      cliqSubFG,
      difflikl.variables,
      difflikl.likelihood;
      graphinit = false,
      tags = [:__LIKELIHOODMESSAGE__; :__UPWARD_DIFFERENTIAL__],
    )
  end

  # listVarByDim = Symbol[]
  # for (label, val) in msgs.belief
  #   push!(listVarByDim, label)
  #   if !exists(tfg, label)
  #     addVariable!(tfg, label, val.variableType)
  #     @debug "New variable added to subfg" _group=:check_addLHDiff #TODO JT remove debug. 
  #   end
  #   initVariable!(tfg, label, manikde!(val))
  # end

  # # list all variables in order of dimension size
  # alreadylist = Symbol[]
  # listDims = getDimension.(getVariable.(tfg,listVarByDim))
  # per = sortperm(listDims, rev=true)
  # listVarDec = listVarByDim[per]
  # listVarAcc = reverse(listVarDec)
  # # add all differential factors (without deconvolution values)
  # for sym1_ in listVarDec
  #   push!(alreadylist, sym1_)
  #   for sym2_ in setdiff(listVarAcc, alreadylist)
  #     nfactype = selectFactorType(tfg, sym1_, sym2_)
  #     # assume default helper function # buildFactorDefault(nfactype)
  #     nfct = nfactype()
  #     afc = addFactor!(tfg, [sym1_;sym2_], nfct, graphinit=false, tags=[:DUMMY;])
  #     # calculate the general deconvolution between variables
  #     pts = solveFactorMeasurements(tfg, afc.label)
  #     newBel = manikde!(getManifold(nfactype), pts[1])
  #     # replace dummy factor with real deconv factor using manikde approx belief measurement
  #     fullFct = nfactype(newBel)
  #     deleteFactor!(tfg, afc.label)
  #     addFactor!( cliqSubFG, [sym1_;sym2_], fullFct, graphinit=false, tags=[:__LIKELIHOODMESSAGE__; :__UPWARD_DIFFERENTIAL__] )
  #   end
  # end

  return tfg
end
# default verbNoun API spec (dest, src)
function addLikelihoodsDifferential!(subfg::AbstractDFG, msgs::LikelihoodMessage)
  return addLikelihoodsDifferential!(msgs, subfg)
end

# child CSM calculates the differential factors that should be sent up
# FIXME, must be renamed and standardized
function addLikelihoodsDifferentialCHILD!(
  cliqSubFG::AbstractDFG,
  seps::Vector{Symbol},
  tfg::AbstractDFG = initfg(
    LocalDFG(; solverParams = SolverParams(; N = getSolverParams(cliqSubFG).N)),
  );
  solveKey::Symbol = :default,
)
  #
  # return list of differential factors the parent should add as part upward partial joint posterior
  retlist = MsgRelativeType()

  # create new local dfg and add all the variables with data
  for label in seps
    if !exists(tfg, label)
      addVariable!(tfg, label, getVariableType(cliqSubFG, label))
      @debug "New variable added to subfg" _group = :check_addLHDiff #TODO JT remove debug. 
    end
    initVariable!(tfg, label, getBelief(cliqSubFG, label, solveKey), solveKey)
  end

  # list all variables in order of dimension size
  alreadylist = Symbol[]
  listDims = getDimension.(getVariable.(tfg, seps))
  per = sortperm(listDims; rev = true)
  listVarDec = seps[per]
  listVarAcc = reverse(listVarDec)
  # add all differential factors (without deconvolution values)
  for sym1_ in listVarDec
    push!(alreadylist, sym1_)
    for sym2_ in setdiff(listVarAcc, alreadylist)
      isHom, ftyps = isPathFactorsHomogeneous(cliqSubFG, sym1_, sym2_)
      # chain of user factors are of the same type
      if isHom
        _sft = selectFactorType(tfg, sym1_, sym2_)
        sft = _sft()
        # only take factors that are homogeneous with the generic relative
        if typeof(sft).name == ftyps[1]
          # assume default helper function # buildFactorDefault(nfactype)
          afc = addFactor!(tfg, [sym1_; sym2_], sft; graphinit = false, tags = [:DUMMY;])
          # calculate the general deconvolution between variables
          pred_X, = approxDeconv(tfg, afc.label, solveKey)  # solveFactorMeasurements
          M = getManifold(_sft)
          e0 = getPointIdentity(M)
          pts = exp.(Ref(M), Ref(e0), pred_X)
          newBel = manikde!(sft, pts)
          # replace dummy factor with real deconv factor using manikde approx belief measurement
          fullFct = _sft(newBel)
          deleteFactor!(tfg, afc.label)
          push!(retlist, (; variables = [sym1_; sym2_], likelihood = fullFct))
        end
      end
    end
  end

  return retlist
end

# use variableList to select a sub-subgraph -- useful for disconnected segments of graph
# NOTE expect msgbeliefs to contain all required keys passed in via special variableList
function _calcCandidatePriorBest(
  subfg::AbstractDFG,
  msgbeliefs::Dict,
  # msgs::LikelihoodMessage,
  variableList::Vector{Symbol} = collect(keys(msgbeliefs)),
)
  #
  ## TODO repackage as new function for wider use
  len = length(variableList)
  dims = Vector{Int}(undef, len)
  syms = Vector{Symbol}(undef, len)
  biAdj = Vector{Int}(undef, len)
  # TODO, not considering existing priors for MsgPrior placement at this time
  # priors = Vector{Int}(undef, len)
  i = 0
  for (label, val) in msgbeliefs
    # skip elements not in variableList
    (label in variableList) ? nothing : continue
    # do calculations based on dimension
    i += 1
    dims[i] = getDimension(val.variableType)
    syms[i] = label
    biAdj[i] = ls(subfg, label) |> length
  end
  # work only with highest dimension variable
  maxDim = maximum(dims)
  dimMask = dims .== maxDim
  mdAdj = biAdj[dimMask]
  pe = sortperm(mdAdj; rev = true) # decending

  # @show variableList, keys(msgbeliefs)
  # @show syms
  return (syms[dimMask])[pe][1]
end

"""
    $SIGNATURES

Generate `MsgPriors` required for upward message joint.  Follows which relative factors ("differentials")
should also be added.

Notes
- Might skip some priors based on `msg.hasPriors`
- This might still be hard to work with, will be clear once engaged in codebase
- TODO obviously much consolidation to do here

Related

[`_findSubgraphsFactorType`](@ref), [`_calcCandidatePriorBest`](@ref)
"""
function _generateSubgraphMsgPriors(
  subfg::AbstractDFG,
  solveKey::Symbol,
  allClasses::Dict{Int, Vector{Symbol}},
  msgbeliefs::Dict,
  msgHasPriors::Bool,
  msgType::MessageType,
)
  #
  priorsJoint = MsgPriorType()

  # 5. find best variable of each of allClasses to place MsgPrior
  for (id, syms) in allClasses
    # if any `jointmsg per variable && !msg.hasPriors`, then dont add a prior
    if (1 == length(syms) || msgHasPriors) && 0 < length(msgbeliefs)
      whichVar = IIF._calcCandidatePriorBest(subfg, msgbeliefs, syms)
      priorsJoint[whichVar] =
        IIF.generateMsgPrior(TreeBelief(getVariable(subfg, whichVar), solveKey), msgType)
    end
  end

  # return the required priors
  return priorsJoint
end

"""
    $SIGNATURES

Generate relative and prior factors that make up the joint msg likelihood.

DevNotes
- Non-standard relative likelihoods will be populated by TAF factors, removing priors assumption.
"""
function _generateMsgJointRelativesPriors(
  cfg::AbstractDFG,
  solveKey::Symbol,
  cliq::TreeClique,
)
  #
  separators = getCliqSeparatorVarIds(cliq)
  jointrelatives = addLikelihoodsDifferentialCHILD!(cfg, separators; solveKey = solveKey)
  allClasses = IIF._findSubgraphsFactorType(cfg, jointrelatives, separators)
  hasPriors = 0 < length(intersect(getCliquePotentials(cliq), lsfPriors(cfg)))

  msgbeliefs = Dict{Symbol, TreeBelief}()
  IIF._buildTreeBeliefDict!(msgbeliefs, cfg, cliq)

  # @show cliq.id, ls(cfg), keys(msgbeliefs), allClasses
  upmsgpriors = IIF._generateSubgraphMsgPriors(
    cfg,
    solveKey,
    allClasses,
    msgbeliefs,
    hasPriors,
    IIF.NonparametricMessage(),
  )

  return _MsgJointLikelihood(; relatives = jointrelatives, priors = upmsgpriors)
end

"""
    $SIGNATURES

Place a single message likelihood prior on the highest dimension variable with highest connectivity in existing subfg.
"""
function addLikelihoodPriorCommon!(
  subfg::AbstractDFG,
  msg::LikelihoodMessage;
  tags::Vector{Symbol} = Symbol[],
)
  #
  tags__ = union(Symbol[:__LIKELIHOODMESSAGE__; :__UPWARD_COMMON__], tags)
  # find if any orphaned variables exist
  for (lbl, msgpr) in msg.jointmsg.priors
    # don't add numerical gauge reference unless absolutely necessary
    if msg.hasPriors || 0 == length(ls(subfg, lbl))
      # finally add the single AbstractPrior from LikelihoodMessage
      addFactor!(subfg, [lbl], msgpr; graphinit = false, tags = tags__)
    end
  end

  # # find max dimension variable, which also has highest biadjacency
  # topCandidate = _calcCandidatePriorBest(subfg, msg.belief)

  # # get prior for top candidate
  # msgPrior = generateMsgPrior(msg.belief[topCandidate], msg.msgType)

  # # get ready
  # tags__ = union(Symbol[:__LIKELIHOODMESSAGE__;:__UPWARD_COMMON__], tags)
  # # finally add the single AbstractPrior from LikelihoodMessage
  # addFactor!(subfg, [topCandidate], msgPrior, graphinit=false, tags=tags__)
end

"""
    $SIGNATURES

Special function to add a few variables and factors to the clique sub graph required for downward solve in CSM.

Dev Notes
- There is still some disparity on whether up and down solves of tree should use exactly the same subgraph...  'between for up and frontal connected for down'
"""
function addDownVariableFactors!(
  dfg::AbstractDFG,
  subfg::InMemoryDFGTypes,
  cliq::TreeClique,
  logger = ConsoleLogger();
  solvable::Int = 1,
)
  #
  # determine which variables and factors needs to be added
  currsyms = ls(subfg)
  allclsyms = getCliqVarsWithFrontalNeighbors(dfg, cliq; solvable = solvable)
  newsyms = setdiff(allclsyms, currsyms)
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.id), newsyms=$newsyms"
  end
  frtls = getCliqFrontalVarIds(cliq)
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.id), frtls=$frtls"
  end
  allnewfcts =
    union(map(x -> findFactorsBetweenFrom(dfg, union(currsyms, newsyms), x), frtls)...)
  newfcts = setdiff(allnewfcts, lsf(subfg))
  with_logger(logger) do
    @info "addDownVariableFactors!, cliq=$(cliq.id), newfcts=$newfcts, allnewfcts=$allnewfcts"
  end

  #TODO solvable?
  DFG.mergeGraph!(subfg, dfg, newsyms, newfcts)

  return newsyms, newfcts
end

"""
    $SIGNATURES

Modify the `subfg::AbstractDFG` to include `msgs` as priors that are used
during clique inference.

Notes
- May be used initialization or inference, in both upward and downward directions.
- `msgs` are identified by variable label `::Symbol`, and my consist of multiple beliefs.
- Message sets from different cliques are identified by clique id `::Int`.
- assume lower limit on number of particles is 5.
- messages from children stored in vector or dict.

DevNotes
- TODO Split dispatch on `dir`, rather than internal `if` statement.

Related

`deleteMsgFactors!`
"""
function addMsgFactors!(
  subfg::AbstractDFG,
  msg::LikelihoodMessage,
  dir::Type{<:MessagePassDirection};
  tags::Vector{Symbol} = Symbol[],
  attemptPriors::Bool = true,
)
  #
  # add messages as priors to this sub factor graph
  msgfcts = DFGFactor[]
  # TODO, expand -- this deconv approach only works for NonparametricMessage at this time.
  if getSolverParams(subfg).useMsgLikelihoods &&
     dir == UpwardPass &&
     msg.msgType isa NonparametricMessage
    #
    if 0 < length(msg.belief)
      # currently only works for nonparametric
      addLikelihoodsDifferential!(subfg, msg)          # :__UPWARD_DIFFERENTIAL__
      if attemptPriors
        # will only be added based on internal tests
        prFcts = addLikelihoodPriorCommon!(subfg, msg)   # :__UPWARD_COMMON__
      end
    end
  else
    svars = DFG.listVariables(subfg)
    tags__ = union(Symbol[:__LIKELIHOODMESSAGE__;], tags)
    dir == DownwardPass ? push!(tags__, :__DOWNWARD_COMMON__) : nothing
    for (msym, belief_) in msg.belief
      if msym in svars
        msgPrior = generateMsgPrior(belief_, msg.msgType)
        fc = addFactor!(subfg, [msym], msgPrior; graphinit = false, tags = tags__)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end

function addMsgFactors_Parametric!(
  subfg::AbstractDFG,
  msg::LikelihoodMessage,
  ::Type{UpwardPass};
  tags::Vector{Symbol} = Symbol[],
  # attemptPriors::Bool = true,
)
  # add differential(relative) message factors

  msgfcts = map(msg.jointmsg.relatives) do difflikl
    addFactor!(
      subfg,
      difflikl.variables,
      difflikl.likelihood;
      graphinit = false,
      tags = union(tags, [:__LIKELIHOODMESSAGE__; :__UPWARD_DIFFERENTIAL__]),
    )
  end

  return msgfcts
end

function addMsgFactors!(
  subfg::AbstractDFG,
  allmsgs::Dict{Int, LikelihoodMessage},
  dir::Type{<:MessagePassDirection};
  tags::Vector{Symbol} = Symbol[],
)
  #
  allfcts = DFGFactor[]
  for (cliqid, msgs) in allmsgs
    # do each dict in array separately
    newfcts = addMsgFactors!(subfg, msgs, dir; tags = tags)
    union!(allfcts, newfcts)
  end
  return allfcts
end

"""
    $SIGNATURES

Delete from the subgraph`::AbstractDFG` prior belief `msgs` that could/would be used
during clique inference.

DevNotes
- TODO make `::Vector{Symbol}` version.
- TODO function taking fcts::Vector{DFGFactor} is unused and replace by the tags version, perhaps we can remove it. 
Related

`addMsgFactors!`
"""
function deleteMsgFactors!(subfg::AbstractDFG, fcts::Vector{DFGFactor})
  #
  for fc in fcts
    deleteFactor!(subfg, fc.label)
  end
end

function deleteMsgFactors!(
  subfg::AbstractDFG,
  tags::Vector{Symbol} = [:__LIKELIHOODMESSAGE__],
)
  # remove msg factors that were added to the subfg
  facs = lsf(subfg; tags = tags)
  deleteFactor!.(subfg, facs)
  return facs
end

## =============================================================================
## Prepare Clique Up or Down Msgs
## =============================================================================

function _buildTreeBeliefDict!(
  msgdict::Dict{Symbol, TreeBelief},
  subfg::AbstractDFG,
  cliq::TreeClique,
  solveKey::Symbol = :default,
  sdims = nothing;    #getCliqVariableMoreInitDims(subfg, cliq, solveKey);
  duplicate::Bool = true,
)
  #
  # TODO better logging
  # with_logger(logger) do
  #   @info "$(now()), prepCliqInitMsgsUp, seps=$seps, sdims=$sdims"
  # end

  seps = getCliqSeparatorVarIds(cliq)
  for vid in seps
    var = DFG.getVariable(subfg, vid)
    var = duplicate ? deepcopy(var) : var
    if isInitialized(var)
      msgdict[var.label] = TreeBelief(var; solvableDim = 1.0) #sdims[var.label])
    end
  end
  return nothing
end

"""
    $SIGNATURES

Prepare the upward inference messages from clique to parent and return as `Dict{Symbol}`.

Notes
- Does not require tree message likelihood factors in subfg.
- Also see #579 regarding elimited likelihoods and priors.

DevNotes
- set `msgs.hasPriors=true` only if a prior occurred here or lower down in tree branch. 
"""
function prepCliqueMsgUp(
  subfg::AbstractDFG,
  cliq::TreeClique,
  solveKey::Symbol,
  status::CliqStatus = getCliqueStatus(cliq);
  logger = ConsoleLogger(),
  duplicate::Bool = true,
  sender = (; id = 0, step = 0),
)
  #
  # get the current clique status
  # sdims = getCliqVariableMoreInitDims(subfg, cliq, solveKey)

  # construct init's up msg to place in parent from initialized separator variables
  hasPriors = 0 < (lsfPriors(subfg) |> length)
  msg = LikelihoodMessage(; sender = sender, status = status, hasPriors = hasPriors)

  _buildTreeBeliefDict!(msg.belief, subfg, cliq, solveKey; duplicate = duplicate)

  # seps = getCliqSeparatorVarIds(cliq)
  # for vid in seps
  #   var = DFG.getVariable(subfg, vid)
  #   var = duplicate ? deepcopy(var) : var
  #   if isInitialized(var)
  #     msg.belief[var.label] = TreeBelief(var, solvableDim=sdims[var.label])
  #   end
  # end

  if getSolverParams(subfg).useMsgLikelihoods
    msg.jointmsg = IIF._generateMsgJointRelativesPriors(subfg, solveKey, cliq)
  end
  # FIXME calculate the new DIFFERENTIAL factors
  # retval = addLikelihoodsDifferentialCHILD!(subfg, getCliqSeparatorVarIds(cliq))
  # msg.jointmsg.relatives = retval

  return msg
end

"""
    $SIGNATURES

Calculate new and then set the down messages for a clique in Bayes (Junction) tree.
"""
function prepCliqueMsgDown(
  subfg::AbstractDFG,
  cliq::TreeClique,
  solveKey::Symbol,
  prntDwnMsgs::LikelihoodMessage,
  logger = ConsoleLogger();
  status::CliqStatus = getCliqueStatus(cliq),
  sender = (; id = cliq.id.value, step = 0),
)
  #
  allvars = getCliqVarIdsAll(cliq)
  allprntkeys = collect(keys(prntDwnMsgs.belief))
  passkeys = intersect(allvars, setdiff(allprntkeys, ls(subfg)))
  remainkeys = setdiff(allvars, passkeys)
  newDwnMsgs = LikelihoodMessage(; sender = sender, status = status)

  # some msgs are just pass through from parent
  for pk in passkeys
    newDwnMsgs.belief[pk] = prntDwnMsgs.belief[pk]
  end

  # set solvable dimensions
  # sdims = getCliqVariableMoreInitDims(subfg, cliq)

  # other messages must be extracted from subfg
  for mk in remainkeys
    setVari = getVariable(subfg, mk)
    if isInitialized(setVari)
      newDwnMsgs.belief[mk] = TreeBelief(setVari, solveKey)  #, solvableDim=sdims[mk] )
    end
  end

  # set the downward keys
  with_logger(logger) do
    @info "cliq $(cliq.id), getSetDownMessagesComplete!, allkeys=$(allvars), passkeys=$(passkeys), msgkeys=$(collect(keys(newDwnMsgs.belief)))"
  end

  return newDwnMsgs
end

## =============================================================================
## Multimessage assemblies from multiple cliques 
## =============================================================================

"""
    $SIGNATURES

Return dictionary of all up belief messages currently in a Bayes `tree`.

"""
function getTreeCliqUpMsgsAll(tree::AbstractBayesTree)
  allUpMsgs = Dict{Int, LikelihoodMessage}()
  for (idx, cliq) in getCliques(tree)
    msgs = getMessageBuffer(cliq).upRx
    merge!(allUpMsgs, msgs)
  end
  return allUpMsgs
end

# TODO @NamedTuple{cliqId::CliqueId{Int}, depth::Int, belief::TreeBelief}
const UpMsgPlotting =
  NamedTuple{(:cliqId, :depth, :belief), Tuple{CliqueId{Int}, Int, TreeBelief}}

"""
    $SIGNATURES

Convert tree up messages dictionary to a new dictionary relative to variables specific messages and their depth in the tree

Notes
- Used in RoMEPlotting
- Return data in `UpMsgPlotting` format.
"""
function stackCliqUpMsgsByVariable(
  tree::AbstractBayesTree,
  tmpmsgs::Dict{Int, LikelihoodMessage},
)
  #
  # start of the return data structure
  stack = Dict{Symbol, Vector{UpMsgPlotting}}()

  # look at all the clique level data
  for (cidx, tmpmsg) in tmpmsgs
    # look at all variables up msg from each clique
    for (sym, belief) in tmpmsg.belief
      # create a new object for a particular variable if hasnt been seen before
      if !haskey(stack, sym)
        # FIXME this is an old message type
        stack[sym] = Vector{UpMsgPlotting}()
      end
      # assemble metadata
      cliq = getClique(tree, cidx)
      #TODO why was the first frontal used? i changed to clique id (unique)
      # frt = getCliqFrontalVarIds(cliq)[1]
      # add this belief msg and meta data to vector of variable entry
      push!(stack[sym], IIF.UpMsgPlotting((cliq.id, getCliqDepth(tree, cliq), belief)))
    end
  end

  return stack
end

"""
    $SIGNATURES

Return dictionary of down messages consisting of all frontal and separator beliefs of this clique.

Notes:
- Fetches numerical results from `subdfg` as dictated in `cliq`.
- return LikelihoodMessage
"""
function getCliqDownMsgsAfterDownSolve(
  subdfg::AbstractDFG,
  cliq::TreeClique,
  solveKey::Symbol;
  status::CliqStatus = NULL,
  sender = (; id = cliq.id.value, step = 0),
)
  #
  # Dict{Symbol, MKD}
  # where the return msgs are contained
  container = LikelihoodMessage(; sender = sender, status = status)

  # go through all msgs one by one
  for sym in getCliqAllVarIds(cliq)
    container.belief[sym] = TreeBelief(getVariable(subdfg, sym), solveKey)
  end

  # return the result
  return container
end
