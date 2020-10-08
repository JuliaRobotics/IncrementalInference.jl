
##==============================================================================
## LEGACY SUPPORT FOR ZMQ IN CAESAR
##==============================================================================

export listSolvekeys

export _evalType

# not sure if and where this is still being used
function _evalType(pt::String)::Type
    @error "_evalType has been deprecated, use DFG serialization methods instead."
    try
        getfield(Main, Symbol(pt))
    catch ex
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
    end
end



##==============================================================================
## Cannot delete until GraphProductOperations usage of this function updated
##==============================================================================

export findRelatedFromPotential

"""
    $(SIGNATURES)

Compute proposal belief on `vertid` through `fct` representing some constraint in factor graph.
Always full dimension variable node -- partial constraints will only influence subset of variable dimensions.
The remaining dimensions will keep pre-existing variable values.

Notes
- fulldim is true when "rank-deficient" -- TODO swap to false (or even float)
"""
function findRelatedFromPotential(dfg::AbstractDFG,
                                  fct::DFGFactor,
                                  varid::Symbol,
                                  N::Int,
                                  dbg::Bool=false;
                                  solveKey::Symbol=:default )
  #
  @warn("findRelatedFromPotential likely to be deprecated, use `lsf` or `productbelief(fg, variableSym, ...) instead`", maxlog=1)

  # assuming it is properly initialized TODO
  ptsbw = evalFactor2(dfg, fct, varid, solveKey=solveKey, N=N, dbg=dbg);
  # determine if evaluation is "dimension-deficient"

  # solvable dimension
  inferdim = getFactorSolvableDim(dfg, fct, varid)
  # zdim = getFactorDim(fct)
  # vdim = getVariableDim(DFG.getVariable(dfg, varid))

  # TODO -- better to upsample before the projection
  Ndim = size(ptsbw,1)
  Npoints = size(ptsbw,2)
  # Assume we only have large particle population sizes, thanks to addNode!
  manis = getManifolds(dfg, varid)
  # manis = getSofttype(DFG.getVariable(dfg, varid)).manifolds # older
  p = AMP.manikde!(ptsbw, manis)
  if Npoints != N # this is where we control the overall particle set size
      p = resample(p,N)
  end
  return (p, inferdim)
end



##==============================================================================
## Delete at end v0.16.x
##==============================================================================

function notifyCliqDownInitStatus!( cliq::TreeClique,
                                    status::Symbol;
                                    logger=ConsoleLogger() )
  #
  @warn("Deprecated, replaced with `prepPutCliqueStatusMsgDwn`")
  cdat = getCliqueData(cliq)
    with_logger(logger) do
    @info "$(now()) $(current_task()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- pre-lock, new $(cdat.initialized)-->$(status)"
  end

  # take lock for atomic transaction
  # lockDwnStatus!(cdat, cliq.index, logger=logger)

  setCliqueStatus!(cdat, status)

  # TODO, should this not send the beliefs aswell??
  msg = LikelihoodMessage(status=status)
  putDwnMsgConsolidated!(cliq, msg)
  # putMsgDwnInitStatus!(cliq, status, logger, msg)
  

  # unlock for others to proceed
  # unlockDwnStatus!(cdat)
  with_logger(logger) do
    @info "$(now()), cliq=$(cliq.index), notifyCliqDownInitStatus! -- unlocked, $(getCliqueStatus(cliq))"
  end

  # flush(logger.stream)

  nothing
end


# function prepareCommonConvWrapper!( ccwl::CommonConvWrapper{F},
#                                     Xi::Vector{DFGVariable},
#                                     solvefor::Symbol,
#                                     N::Int;
#                                     solveKey::Symbol=:default  ) where {F <: AbstractRelative}
#   #
#   prepareCommonConvWrapper!(F, ccwl, Xi, solvefor, N, solveKey=solveKey)
# end

# function computeAcrossHypothesis!(ccwl::Union{CommonConvWrapper{F},CommonConvWrapper{Mixture{N_,F,S,T}}},
#                                   allelements,
#                                   activehypo,
#                                   certainidx::Vector{Int},
#                                   sfidx::Int,
#                                   maxlen::Int,
#                                   maniAddOps::Tuple;
#                                   spreadNH::Real=3.0  ) where {N_,F<:AbstractRelative,S,T}
#   #
#   computeAcrossHypothesis!(F,ccwl,allelements,activehypo,certainidx,sfidx, maxlen,maniAddOps,spreadNH=spreadNH)
# end


# function computeAcrossHypothesis!(ccwl::CommonConvWrapper{Mixture{N_,F,S,T}},
#                                   allelements,
#                                   activehypo,
#                                   certainidx::Vector{Int},
#                                   sfidx::Int,
#                                   maxlen::Int,
#                                   maniAddOps::Tuple;
#                                   spreadNH::Real=3.0  ) where {N_,F<:AbstractRelative,S,T}
#   #
#   computeAcrossHypothesis!(F,ccwl,allelements,activehypo,certainidx,sfidx, maxlen,maniAddOps,spreadNH=spreadNH)
# end


@deprecate numericRootGenericRandomizedFnc!(w...;kw...) numericSolutionCCW!(w...;kw...)


# function numericRootGenericRandomizedFnc!(ccwl::CommonConvWrapper{Mixture{N,F,S,T}};
#                                           perturb::Float64=1e-10,
#                                           testshuffle::Bool=false ) where 
#                                               {N,F<:AbstractRelative,S,T <: Tuple}
#   #
#   _numericSolutionCCW!(F, ccwl,perturb=perturb, testshuffle=testshuffle)
# end


# function numericRootGenericRandomizedFnc!(ccwl::CommonConvWrapper{F};
#                                           perturb::Float64=1e-10,
#                                           testshuffle::Bool=false ) where 
#                                               {F <: AbstractRelative}
#   #
#   _numericSolutionCCW!(F, ccwl, perturb=perturb, testshuffle=testshuffle)
# end


@deprecate MixtureRelative(w...; kw...) Mixture(w...; kw...)

@deprecate MixturePrior(w...; kw...) Mixture(Prior, w...; kw...)


# function areSiblingsRemaingNeedDownOnly(tree::AbstractBayesTree,
#   cliq::TreeClique  )::Bool
#   #
#   stillbusylist = [:null; :initialized;]
#   prnt = getParent(tree, cliq)
#   if length(prnt) > 0
#     for si in getChildren(tree, prnt[1])
#       # are any of the other siblings still busy?
#       if si.index != cliq.index && getCliqueStatus(si) in stillbusylist
#         noOneElse = false
#       end
#     end
#   end
  
#   # nope, everybody is waiting for something to change -- proceed with forcing a cliq solve
#   return true
# end


# """
#     $SIGNATURES

# Return true if both, i.) this clique requires more downward information, ii.) more
# downward message information could potentially become available.

# Notes
# - Delay initialization to the last possible moment.

# Dev Notes:

# Determine clique truely isn't able to proceed any further:
# - should be as self reliant as possible (using clique's status as indicator)
# - change status to :mustinitdown if have only partial beliefs so far:
#   - combination of status, while partials belief siblings are not :mustinitdown
# """
# function getCliqSiblingsPartialNeeds(tree::AbstractBayesTree,
#                                       cliq::TreeClique,
#                                       #  prnt,
#                                       dwinmsgs::LikelihoodMessage;
#                                       logger=ConsoleLogger())
#   #
#   # which incoming messages are partials
#   hasPartials = Dict{Symbol, Int}()
#   for (sym, tmsg) in dwinmsgs.belief
#     # assuming any down message per label that is not partial voids further partial consideration
#     if sum(tmsg.inferdim) > 0
#       if !haskey(hasPartials, sym)
#         hasPartials[sym] = 0
#       end
#       hasPartials[sym] += 1
#     end
#   end
#   partialKeys = collect(keys(hasPartials))
  
#   ## determine who might be able to help init this cliq
#   # check sibling separator sets against this clique's separator
#   sibs = getCliqSiblings(tree, cliq)
  
#   with_logger(logger) do
#     @info "getCliqSiblingsPartialNeeds -- CHECK PARTIAL"
#   end
#   # identify which cliques might have useful information
#   localsep = getCliqSeparatorVarIds(cliq)
#   seps = Dict{Int, Vector{Symbol}}()
#   for si in sibs
#     # @show getLabel(si)
#     mighthave = intersect(getCliqSeparatorVarIds(si), localsep)
#     if length(mighthave) > 0
#       seps[si.index] = mighthave
#       if getCliqueStatus(si) in [:initialized; :null; :needdownmsg]
#         # partials treated special -- this is slightly hacky
#         if length(intersect(localsep, partialKeys)) > 0 && length(mighthave) > 0
#           # this sibling might have info to delay about
#           setCliqDrawColor(cliq,"magenta")
#           return true
#         end
#       end
#     end
#   end
#   # determine if those cliques will / or will not be able to provide more info
#   # when does clique change to :mustinitdown
#   # default
#   return false
# end


# """
#     $SIGNATURES

# Return true if this clique's down init should be delayed on account of prioritization among sibling separators.

# Notes
# - process described in issue #344

# Dev Notes
# - not priorizing order yet (TODO), just avoiding unsolvables at this time.
# - Very closely related to getCliqSiblingsPartialNeeds -- refactor likely (NOTE).
# - should precompute `allinters`.
# """
# function getSiblingsDelayOrder(tree::AbstractBayesTree,
#                                 cliq::TreeClique,
#                                 #  prnt,
#                                 dwnkeys::Vector{Symbol}; # dwinmsgs::LikelihoodMessage;
#                                 logger=ConsoleLogger())
#   # when is a cliq upsolved
#   solvedstats = Symbol[:upsolved; :marginalized; :uprecycled]

#   # safety net double check
#   cliqst = getCliqueStatus(cliq)
#   if cliqst in solvedstats
#     with_logger(logger) do
#       @warn "getSiblingsDelayOrder -- clique status should not be here with a solved cliqst=$cliqst"
#     end
#     return false
#   end

#   # get siblings separators
#   sibs = getCliqSiblings(tree, cliq, true)
#   ids = map(s->s.index, sibs)
#   len = length(sibs)
#   sibidx = collect(1:len)[ids .== cliq.index][1]
#   seps = getCliqSeparatorVarIds.(sibs)
#   lielbls = setdiff(ids, cliq.index)
#   # get intersect matrix of siblings (should be exactly the same across siblings' csm)
#   allinters = Array{Int,2}(undef, len, len)
#   dwninters = Vector{Int}(undef, len)
#   with_logger(logger) do
#     @info "getSiblingsDelayOrder -- number siblings=$(len), sibidx=$sibidx"
#   end

#   # sum matrix with all "up solved" rows and columns eliminated
#   fill!(allinters, 0)
#   for i in 1:len
#     for j in i:len
#       if i != j
#         allinters[i,j] = length(intersect(seps[i],seps[j]))
#       end
#     end
#     dwninters[i] = length(intersect(seps[i], dwnkeys))
#   end

#   # sum "across/over" rows, then columns (i.e. visa versa "along" columns then rows)
#   rows = sum(allinters, dims=1)
#   cols = sum(allinters, dims=2)

#   with_logger(logger) do
#       @info "getSiblingsDelayOrder -- allinters=$(allinters)"
#       @info "getSiblingsDelayOrder -- rows=$(rows)"
#       @info "getSiblingsDelayOrder -- rows=$(cols)"
#   end

#   # is this clique a non-zero row -- i.e. sum across columns? if not, no further special care needed
#   if cols[sibidx] == 0
#     with_logger(logger) do
#       @info "getSiblingsDelayOrder -- cols[sibidx=$(sibidx))] == 0, no special care needed"
#     end
#     return false
#   end

#   # now determine if initializing from below or needdownmsg
#   if cliqst in Symbol[:needdownmsg;]
#     # be super careful about delay (true) vs pass (false) at this point -- might be partial too TODO
#     # return true if delay beneficial to initialization accuracy

#     # find which siblings this cliq epends on
#     symm = allinters + allinters'
#     maskcol = symm[:,sibidx] .> 0
#     # lenm = length(maskcol)
#     stat = Vector{Symbol}(undef, len)
#     stillbusymask = fill(false, len)

#     flush(logger.stream)

#     # get each sibling status (entering atomic computation segment -- until wait command)
#     stat .= getCliqueStatus.(sibs) #[maskcol]

#     ## (long down chain case)
#     # need different behaviour when all remaining siblings are blocking with :needdownmsg
#     remainingmask = stat .== :needdownmsg
#     if sum(remainingmask) == length(stat)
#       with_logger(logger) do
#         @info "getSiblingsDelayOrder -- all blocking: sum(remainingmask) == length(stat), stat=$stat"
#       end
#       # pick sibling with most overlap in down msgs from parent
#       # list of similar length siblings
#       candidates = dwninters .== maximum(dwninters)
#       if candidates[sibidx]
#         # must also pick minimized intersect with other remaing siblings
#         maxcan = collect(1:len)[candidates]
#         with_logger(logger) do
#           @info "getSiblingsDelayOrder -- candidates=$candidates, maxcan=$maxcan, rows=$rows"
#         end
#         if rows[sibidx] == minimum(rows[maxcan])
#           with_logger(logger) do
#             @info "getSiblingsDelayOrder -- FORCE DOWN INIT SOLVE ON THIS CLIQUE: $(cliq.index), $(getLabel(cliq))"
#           end
#           return false
#         end
#       end
#       with_logger(logger) do
#         @info "getSiblingsDelayOrder -- not a max and should block"
#       end
#       return true
#     end

#     # still busy solving on branches, so potential to delay
#     for i in 1:len
#       stillbusymask[i] = maskcol[i] && !(stat[i] in solvedstats)
#     end
#     with_logger(logger) do
#         @info "getSiblingsDelayOrder -- busy solving:"
#         @info "maskcol=$maskcol"
#         @info "stillbusy=$stillbusymask"
#     end

#     # Too blunt -- should already have returned false by this point perhaps
#     if sum(stillbusymask) > 0
#       # yes something to delay about
#       with_logger(logger) do
#         @info "getSiblingsDelayOrder -- yes delay,"
#         @info "stat=$stat"
#         @info "symm=$symm"
#       end
#       return true
#     end
#   end

#   with_logger(logger) do
#     @info "getSiblingsDelayOrder -- default will not delay"
#   end
#   flush(logger.stream)
#   # carry over default from partial init process
#   return false
# end


@deprecate getMsgsUpChildren(x...;kw...) fetchMsgsUpChildren(x...;kw...)

@deprecate getMsgUpThis(x...;kw...) fetchMsgUpThis(x...;kw...)

@deprecate printCliqHistorySequential(x...;kw...) printCSMHistorySequential(x...;kw...)

# """
#     $SIGNATURES

# WIP #459 dwnMsg consolidation towards blocking cliq that `:needdwninit` to wait on parent `:initialized` dwn message.

# Notes
# - State machine function nr.6e

# DevNotes
# - Seems really unnecessary
# - Separated out during #459 dwnMsg consolidation
# - Should only happen in long downinit chains below parent that needed dwninit
# - TODO figure out whats different between this and 8c
# """
# function slowOnPrntAsChildrNeedDwn_StateMachine(csmc::CliqStateMachineContainer)
#   # do actual fetch
#   prtmsg = fetchDwnMsgConsolidated(getParent(csmc.tree, csmc.cliq)[1]).status

#   # FIXME WHY THIS???
#   # go to 7
#   return determineCliqNeedDownMsg_StateMachine
# end


# """
# $SIGNATURES

# WIP to resolve 459 dwnMsg consolidation.  This is partly doing some kind of downsolve but seems out of place.

# Notes
# - State machine function nr. 10a

# DevNotes
# - FIXME, resolve/consolidate with 8c?
# """
# function wipRedirect459Dwn_StateMachine(csmc::CliqStateMachineContainer)
#   infocsm(csmc, "10a, canCliqDownSolve_StateMachine, going to block on parent.")
#   prnt = getParent(csmc.tree, csmc.cliq)

#   # block here until parent is downsolved
#   setCliqDrawColor(csmc.cliq, "turquoise")
#   # this part is a pull model #674
#   # while 
#   prntst = fetchDwnMsgConsolidated(prnt[1]).status
#   if prntst != :downsolved
#     wait(getSolveCondition(prnt[1]))
#   end
#   # blockMsgDwnUntilStatus(prnt[1], :downsolved)
#   # blockCliqUntilParentDownSolved(, logger=csmc.logger)

#   # yes, continue with downsolve
#   # prntst = getCliqueStatus(prnt[1])
#   infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, parent status=$prntst.")
#   if prntst != :downsolved
#     infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, going around again.")
#     return canCliqDownSolve_StateMachine
#   end

#   infocsm(csmc, "10a, wipRedirect459Dwn_StateMachine, going for down solve.")
#   # go to 11
#   return doCliqDownSolve_StateMachine
# end



"""
$(TYPEDEF)

TO BE DEPRECATED AND CONSOLIDATED
"""
mutable struct DownReturnBPType
  dwnMsg::LikelihoodMessage
  dbgDwn::DebugCliqMCMC
  IDvals::Dict{Symbol,TreeBelief}
  keepdwnmsgs::LikelihoodMessage
end


#NOTE select type for development
# emptyBayesTree() = BayesTree()
# emptyBayesTree() = MetaBayesTree()

function dwnPrepOutMsg( fg::AbstractDFG,
                        cliq::TreeClique,
                        dwnMsgs::Array{LikelihoodMessage,1},
                        d::Dict{Symbol, T},
                        logger=ConsoleLogger()) where T
  @warn("dwnPrepOutMsg is deprecated")                      
  # pack all downcoming conditionals in a dictionary too.
  with_logger(logger) do
    if cliq.index != 1 #TODO there may be more than one root
      @info "Dwn msg keys $(keys(dwnMsgs[1].belief))"
      @info "fg vars $(ls(fg))"
    end # ignore root, now incoming dwn msg
  end
  m = LikelihoodMessage()
  i = 0
  for vid in getCliqueData(cliq).frontalIDs
    m.belief[vid] = deepcopy(d[vid]) # TODO -- not sure if deepcopy is required
  end
  for cvid in getCliqueData(cliq).separatorIDs
    i+=1
    # TODO -- convert to points only since kde replace by rkhs in future
    m.belief[cvid] = deepcopy(dwnMsgs[1].belief[cvid]) # TODO -- maybe this can just be a union(,)
  end
  return m
end



"""
    $SIGNATURES

Perform Chapman-Kolmogorov transit integral approximation for `cliq` in downward pass direction.

Notes
- Only update frontal variables of the clique.
"""
function downGibbsCliqueDensity(fg::AbstractDFG,
                                cliq::TreeClique,
                                dwnMsgs::Array{LikelihoodMessage,1},
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger() )
  
  @warn("downGibbsCliqueDensity(fg, cliq, dwnMsgs::Array{LikelihoodMessage,1}, N, MCMCIter, dbg, usemsgpriors, logger) is deprecated")                                   
  # TODO standardize function call to have similar stride to upGibbsCliqueDensity
  # @info "down"
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- going for down fmcmc"
  end
  fmmsgs = usemsgpriors ? Array{LikelihoodMessage,1}() : dwnMsgs
  frtls = getFrontals(cliq)

  # TODO, do better check if there is structure between multiple frontals
  niters = length(frtls) == 1 ? 1 : MCMCIter
  # TODO standize with upsolve and variable solver order
  mcmcdbg, d = fmcmc!(fg, cliq, fmmsgs, frtls, N, niters, dbg)
  m = dwnPrepOutMsg(fg, cliq, dwnMsgs, d, logger)

  outmsglbl = Dict{Symbol, Int}()
  if dbg
    for (ke, va) in m.belief
      outmsglbl[Symbol(fg.g.vertices[ke].label)] = ke
    end
  end

  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- convert to BallTreeDensities."
  end

  # Always keep dwn messages in cliq data
  dwnkeepmsgs = LikelihoodMessage()
  for (msgsym, val) in m.belief
    dwnkeepmsgs.belief[msgsym] = convert(Tuple{BallTreeDensity,Float64}, val)
  end
  #FIXME this does not exist so looks like this function is never called!
  setDwnMsg!(cliq, dwnkeepmsgs)

  # down solving complete, set flag
  getCliqueData(cliq).downsolved = true

  mdbg = !dbg ? DebugCliqMCMC() : DebugCliqMCMC(mcmcdbg, m, outmsglbl, CliqGibbsMC[])
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- finished."
  end
  return DownReturnBPType(m, mdbg, d, dwnkeepmsgs)
end


function downGibbsCliqueDensity(fg::AbstractDFG,
                                cliq::TreeClique,
                                dwnMsgs::LikelihoodMessage,
                                N::Int=100,
                                MCMCIter::Int=3,
                                dbg::Bool=false,
                                usemsgpriors::Bool=false,
                                logger=ConsoleLogger() )
  #
  @warn("downGibbsCliqueDensity is deprecated")
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- convert BallTreeDensities to LikelihoodMessage."
  end
  ind = Dict{Symbol, TreeBelief}()
  sflbls = listVariables(fg)
  for (lbl, bel) in dwnMsgs.belief
    if lbl in sflbls
      ind[lbl] = TreeBelief(bel[1], bel[2], getSofttype(getVariable(fg, lbl)))
    end
  end
  ndms = LikelihoodMessage[LikelihoodMessage(ind);]
  with_logger(logger) do
    @info "cliq=$(cliq.index), downGibbsCliqueDensity -- call with LikelihoodMessage."
  end
  downGibbsCliqueDensity(fg, cliq, ndms, N, MCMCIter, dbg, usemsgpriors, logger)
end


## TODO IS THIS DEPRECATED - never called 
"""
    $(SIGNATURES)

Update cliq `cliqID` in Bayes (Juction) tree `bt` according to contents of `ddt` -- intended use is to update main clique after a downward belief propagation computation has been completed per clique.
"""
function updateFGBT!( fg::AbstractDFG,
                      bt::AbstractBayesTree,
                      cliqID::Int,
                      drt::DownReturnBPType;
                      dbg::Bool=false,
                      fillcolor::String="",
                      logger=ConsoleLogger()  )
    
    @warn("updateFGBT!(fg, bt, cliqID, drt::DownReturnBPType) is deprecated")
    cliq = getClique(bt, cliqID)
    # if dbg
    #   cliq.attributes["debugDwn"] = deepcopy(drt.dbgDwn)
    # end
    #FIXME this does not exist so looks like this function is never called!
    setDwnMsg!(cliq, drt.keepdwnmsgs)
    # TODO move to drawTree
    if fillcolor != ""
      setCliqDrawColor(cliq, fillcolor)
    end
    with_logger(logger) do
      for (sym, emsg) in drt.IDvals
        #TODO -- should become an update call
        updvert = DFG.getVariable(fg, sym)
        # TODO -- not sure if deepcopy is required , updatePPE=true)
        @info "updateFGBT, DownReturnBPType, sym=$sym, current inferdim val=$(getVariableInferredDim(updvert))"
        setValKDE!(updvert, deepcopy(emsg) )
        updvert = DFG.getVariable(fg, sym)
        @info "updateFGBT, DownReturnBPType, sym=$sym, inferdim=$(emsg.inferdim), newval=$(getVariableInferredDim(updvert))"
      end
    end
    nothing
end


@deprecate prepCliqInitMsgsUp(x...;kw...) prepCliqueMsgUpConsolidated(x...;kw...)
@deprecate getSetDownMessagesComplete!(x...;kw...) prepSetCliqueMsgDownConsolidated!(x...;kw...)

@deprecate getMsgDwnChannel(tree::AbstractBayesTree, edge) getDwnMsgConsolidated(tree, edge)

# export MixtureLinearConditional

@deprecate MixtureLinearConditional(Z::AbstractVector{<:SamplableBelief}, C::DiscreteNonParametric)  Mixture(LinearRelative, Z, C)


"""
    $(SIGNATURES)

Add all potentials associated with this variable in `cliq` to `dens`.
"""
function packFromLocalPotentials!(dfg::AbstractDFG,
                                  dens::Vector{BallTreeDensity},
                                  wfac::Vector{Symbol},
                                  cliq::TreeClique,
                                  vsym::Symbol,
                                  N::Int,
                                  dbg::Bool=false )
  #
  @warn("packFromLocalPotentials! is obsolete, use `productbelief(fg, variableSym, :)`")

  inferdim = 0.0
  for idfct in getCliqueData(cliq).potentials
    !(exists(dfg, idfct)) && (@warn "$idfct not in clique $(cliq.index)" continue)
    fct = DFG.getFactor(dfg, idfct)
    data = getSolverData(fct)
    # skip partials here, will be caught in packFromLocalPartials!
    if length( findall(getVariableOrder(fct) .== vsym) ) >= 1 && !data.fnc.partial
    # if length( findall(data.fncargvID .== vsym) ) >= 1 && !data.fnc.partial
      p, isinferdim = findRelatedFromPotential(dfg, fct, vsym, N, dbg )
      push!(dens, p)
      push!(wfac, fct.label)
      inferdim += isinferdim
    end
  end

  # return true if at least one of the densities was full dimensional (used for tree based initialization logic)
  return inferdim::Float64
end


function packFromLocalPartials!(fgl::AbstractDFG,
                                partials::Dict{Int, Vector{BallTreeDensity}},
                                cliq::TreeClique,
                                vsym::Symbol,
                                N::Int,
                                dbg::Bool=false  )
  #
  @warn("packFromLocalPartials! is obsolete, use `productbelief(fg, variableSym, :)`")

  for idfct in getCliqueData(cliq).potentials
    !(exists(fgl, idfct)) && (@warn "$idfct not in clique $(cliq.index)" continue)
    vert = DFG.getFactor(fgl, idfct)
    data = getSolverData(vert)
    if length( findall(getVariableOrder(vert) .== vsym) ) >= 1 && data.fnc.partial
      p, = findRelatedFromPotential(fgl, vert, vsym, N, dbg)
      pardims = data.fnc.usrfnc!.partial
      for dimnum in pardims
        if haskey(partials, dimnum)
          push!(partials[dimnum], marginal(p,[dimnum]))
        else
          partials[dimnum] = BallTreeDensity[marginal(p,[dimnum])]
        end
      end
    end
  end
  nothing
end


"""
    $SIGNATURES

Previous generation function converting LikelihoodMessage into Vector of densities before inference product.

DevNotes
- #913 and consolidation with CSM and `csmc.cliqSubFg`
- FIXME FIXME make sure this is using the proper densities/potentials/likelihoods/factors stored in `csmc.cliqSubFg`.
"""
function packFromIncomingDensities!(dens::Vector{BallTreeDensity},
                                    wfac::Vector{Symbol},
                                    vsym::Symbol,
                                    inmsgs::Array{LikelihoodMessage,1},
                                    manis::T )::Float64 where {T <: Tuple}
  #
  @warn("packFromIncomingDensities! is obsolete, use `productbelief(fg, variableSym, :)`")
  inferdim = 0.0
  for m in inmsgs
    for psym in keys(m.belief)
      if psym == vsym
        pdi = m.belief[vsym]
        push!(dens, manikde!(pdi.val, pdi.bw[:,1], getManifolds(pdi)) )
        push!(wfac, :msg)
        inferdim += pdi.inferdim
      end
      # TODO -- we can inprove speed of search for inner loop
    end
  end

  # return true if at least one of the densities was full dimensional (used for tree based initialization logic)
  return inferdim
end


"""
    $(SIGNATURES)

Perform one step of the minibatch clique Gibbs operation for solving the Chapman-Kolmogov
trasit integral -- here involving separate approximate functional convolution and
product operations.
"""
function cliqGibbs( dfg::AbstractDFG,
                    cliq::TreeClique,
                    vsym::Symbol,
                    inmsgs::Array{LikelihoodMessage,1},
                    N::Int,
                    dbg::Bool,
                    manis::Tuple,
                    logger=ConsoleLogger()  )
  #
  @warn("cliqGibbs is deprecated, use productbelief(fg, variableSym, :) instead")
  # several optimizations can be performed in this function TODO
  
  inferdim = 0.0
  # consolidate NBPMessages and potentials
  dens = Array{BallTreeDensity,1}()
  partials = Dict{Int, Vector{BallTreeDensity}}()
  wfac = Vector{Symbol}()
  # figure out which densities to use
  inferdim += packFromIncomingDensities!(dens, wfac, vsym, inmsgs, manis)
  inferdim += packFromLocalPotentials!(dfg, dens, wfac, cliq, vsym, N)
  packFromLocalPartials!(dfg, partials, cliq, vsym, N, dbg)

  potprod = !dbg ? nothing : PotProd(vsym, getVal(dfg,vsym), Array{Float64,2}(undef, 0,0), dens, wfac)
      # pts,inferdim = predictbelief(dfg, vsym, useinitfct)  # for reference only
  pGM = productbelief(dfg, vsym, dens, partials, N, dbg=dbg, logger=logger )

  if dbg  potprod.product = pGM  end

  # @info " "
  return pGM, potprod, inferdim
end



export FullExploreTreeType, ExploreTreeType

"""
$(TYPEDEF)
"""
mutable struct FullExploreTreeType{T, T2, T3 <:InMemoryDFGTypes}
  fg::T3
  bt::T2
  cliq::TreeClique
  prnt::T
  sendmsgs::Vector{LikelihoodMessage}
end

const ExploreTreeType{T} = FullExploreTreeType{T, BayesTree}
const ExploreTreeTypeLight{T} = FullExploreTreeType{T, Nothing}


function ExploreTreeType( fgl::G,
                          btl::AbstractBayesTree,
                          vertl::TreeClique,
                          prt::T,
                          msgs::Array{LikelihoodMessage,1} ) where {G <: AbstractDFG, T}
  #
  ExploreTreeType{T}(fgl, btl, vertl, prt, msgs)
end

@deprecate upGibbsCliqueDensity(inp::FullExploreTreeType{T,T2},
                                N::Int=100,
                                dbg::Bool=false,
                                iters::Int=3,
                                logger=ConsoleLogger()  ) where {T, T2} upGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, dbg, iters, logger)


export initVariable!

"""
    $(SIGNATURES)

Initialize the belief of a variable node in the factor graph struct.
"""
function initVariable!(fgl::AbstractDFG,
                        sym::Symbol;
                        N::Int=100  )
  #
  @warn "initVariable! has been displaced by doautoinit! or initManual! -- might be revived in the future"

  vert = getVariable(fgl, sym)
  belief,b,c,d,infdim  = localProduct(fgl, sym, N=N)
  setValKDE!(vert, belief)

  nothing
end


#global pidx
global pidx = 1
global pidl = 1
global pidA = 1
global thxl = nprocs() > 4 ? floor(Int,nprocs()*0.333333333) : 1

# upploc to control processes done local to this machine and separated from other
# highly loaded processes. upploc() should be used for dispatching short and burst
# of high bottle neck computations. Use upp2() for general multiple dispatch.
function upploc()
    @warn "upploc is deprecated, use WorkerPool instead"
    global pidl, thxl
    N = nprocs()
    pidl = (thxl > 1 ? ((pidl < thxl && pidl != 0 && thxl > 2) ? (pidl+1)%(thxl+1) : 2) : (N == 1 ? 1 : 2))
    return pidl
end

# upp2() may refer to processes on a different machine. Do not associate upploc()
# with processes on a separate computer -- this will be become more complicated when
# remote processes desire their own short fast 'local' burst computations. We
# don't want all the data traveling back and forth for shorter jobs
function upp2()
  @warn "upploc is deprecated, use WorkerPool instead"
  global pidx, thxl
  N = nprocs()
  pidx = (N > 1 ? ((pidx < N && pidx != 0) ? (pidx+1)%(N+1) : thxl+1) : 1) #2 -- thxl+1
  return pidx
end

function uppA()
  @warn "upploc is deprecated, use WorkerPool instead"
  global pidA
  N = nprocs()
  pidA = (N > 1 ? ((pidA < N && pidA != 0) ? (pidA+1)%(N+1) : 2) : 1) #2 -- thxl+1
  return pidA
end




@deprecate wipeBuildNewTree!(dfg::AbstractDFG; kwargs...) resetBuildTree!(dfg; kwargs...)

# getSample(s::MixtureRelative, N::Int=1) = (reshape.(rand.(s.Z, N),1,:)..., rand(s.C, N))

# @deprecate (MixturePrior{T}(z::NTuple{N,<:SamplableBelief}, c::Union{<:Distributions.DiscreteNonParametric, NTuple{N,<:Real}, <:AbstractVector{<:Real}} ) where {T,N}) MixturePrior(z,c)

@deprecate LinearConditional(N::Int=1) LinearRelative{N}(LinearAlgebra.I)
# @deprecate LinearConditional(x::SamplableBelief) LinearRelative(x)
@deprecate LinearConditional(x...) LinearRelative(x...)

# function LinearConditional{N, T}(x...) where {N,T}
#   @warn("LinearConditional{N, T} is deprecated, use LinearRelative instead")
#   LinearRelative{N,T}(x...)
# end

@deprecate PackedLinearConditional(x...) PackedLinearRelative(x...)



# =============================================================================
# Clique condition locks

@deprecate lockUpStatus!(x...) ()->nothing

@deprecate unlockUpStatus!(cdat::BayesTreeNodeData) ()->nothing
unlockUpStatus!(cliq::TreeClique) = unlockUpStatus!(getCliqueData(cliq))

@deprecate lockDwnStatus!(cdat::BayesTreeNodeData, idx::Int=1; logger=ConsoleLogger()) ()->nothing

@deprecate unlockDwnStatus!(cdat::BayesTreeNodeData) ()->nothing




function iifdepwarn(msg, funcsym; maxlog=nothing)
  @logmsg(
      Base.CoreLogging.Warn,
      msg,
      _module=begin
          bt = backtrace()
          frame, caller = Base.firstcaller(bt, funcsym)
          # TODO: Is it reasonable to attribute callers without linfo to Core?
          caller.linfo isa Core.MethodInstance ? caller.linfo.def.module : Core
      end,
      _file=String(caller.file),
      _line=caller.line,
      _id=(frame,funcsym),
      _group=:iifdepwarn,
      caller=caller,
      short_stacktrace=stacktrace(bt)[7:9],
      maxlog=maxlog
  )
  nothing
end

function Base.getproperty(obj::BayesTreeNodeData, sym::Symbol)
  if sym == :dwnMsg
    # iifdepwarn("#459 get dwnMsg", :getproperty)
  elseif sym == :downInitMsg
    # iifdepwarn("#459 get downInitMsg", :getproperty)
  elseif sym == :initDownChannel
    # iifdepwarn("#459 get initDownChannel", :getproperty)
  end
  return getfield(obj, sym)
end

function Base.setproperty!(obj::BayesTreeNodeData, sym::Symbol, val)
  if sym == :dwnMsg
    # iifdepwarn("#459 set dwnMsg", :setproperty!)
  elseif sym == :downInitMsg
    # iifdepwarn("#459 set downInitMsg", :setproperty!)
  elseif sym == :initDownChannel
    # iifdepwarn("#459 set initDownChannel", :setproperty!)
  end
  return setfield!(obj, sym, convert(fieldtype(typeof(obj), sym), val))
end


# ============================================================================
# .initDownChannel, MUST BE REMOVED

# ============================================================================
# .downInitMsg, MUST BE REMOVED


@deprecate putMsgDwnInitChannel!(btnd::BayesTreeNodeData, msg::LikelihoodMessage) putDwnMsgConsolidated!(btnd, msg)
@deprecate getMsgDwnInitChannel_(btnd::BayesTreeNodeData) getDwnMsgConsolidated(btnd)
# getMsgDwnInitChannel_(btnd::BayesTreeNodeData) = btnd.initDownChannel

function getMsgDwnInitChannel_(cliq::TreeClique)
  @warn("getMsgDwnInitChannel_ is deprecated, use getDwnMsgConsolidated instead.", maxlog=1)
  getMsgDwnInitChannel_(getCliqueData(cliq))
end
export fetchMsgDwnInit
fetchMsgDwnInit(cliq::TreeClique) = fetch(getMsgDwnInitChannel_(cliq))


# FIXME OLD must be consolidated as part of 459
function putMsgDwnInitStatus!(cliq::TreeClique, status::CliqStatus, logger=ConsoleLogger(), msg=LikelihoodMessage(status=status))
  @warn("putMsgDwnInitStatus! is deprecated, use putDwnMsgConsolidated! instead")
  cdat = getCliqueData(cliq)
  cdc = getDwnMsgConsolidated(cdat)
  # cdc = getMsgDwnInitChannel_(cdat)
    if isready(cdc)
      content = take!(cdc)
      with_logger(logger) do
        @info "dumping stale cliq=$(cliq.index) status message $(content), replacing with $(status)"
      end
    end
  put!(cdc, msg)
  notify(getSolveCondition(cliq))
    # FIXME hack to mitigate old race condition
    sleep(0.1)
    notify(getSolveCondition(cliq))

  nothing
end

function getfetchCliqueInitMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing)
  @debug "getfetchCliqueInitMsgDown from=$(from)"
  return cdata.downInitMsg
end
# getMsgDwnThisInit(cliq::TreeClique) = getMsgDwnThisInit(getCliqueData(cliq)) # WHAT ???

function putCliqueInitMsgDown!(cdata::BayesTreeNodeData, initmsg::LikelihoodMessage)
  @warn("putCliqueInitMsgDown! is deprecated, use putDwnMsgConsolidated instead")
  cdata.downInitMsg = initmsg
  nothing
end


@deprecate fetchMsgDwnThis(cliq::TreeClique) fetchDwnMsgConsolidated(cliq)

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
fetchMsgDwnThis(csmc::CliqStateMachineContainer) = fetchMsgDwnThis(csmc.cliq)
fetchMsgDwnThis(btl::AbstractBayesTree, sym::Symbol) = fetchMsgDwnThis(getClique(btl, sym))
# function fetchMsgDwnThis(cliql::TreeClique)
#   @error("fetchMsgDwnThis is deprecated, use getDwnMsgConsolidated instead")
#  # fetchDwnMsgConsolidated(cliql)
  # getCliqueData(cliql).dwnMsg
# end


@deprecate putMsgDwnThis!(x...; kw...) putDwnMsgConsolidated!(x...; kw...)

"""
$(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""  
putMsgDwnThis!(csmc::CliqStateMachineContainer, msgs::LikelihoodMessage) = putMsgDwnThis!(csmc.cliq, msgs)  # NOTE, old, csmc.msgsDown = msgs
# function putMsgDwnThis!(cdata::BayesTreeNodeData, msg::LikelihoodMessage; from::Symbol=:nothing)
#   @error("putMsgDwnThis! is deprecated, use putDwnMsgConsolidated instead")
#   @debug "putMsgDwnThis! from=$(from)"
#   cdata.dwnMsg = msg
# end  
# function putMsgDwnThis!(cliql::TreeClique, msgs::LikelihoodMessage)
#   getCliqueData(cliql).dwnMsg = msgs
# end  
  




"""
    $SIGNATURES

Do down initialization calculations, loosely translates to solving Chapman-Kolmogorov
transit integral in downward direction.

Notes
- State machine function nr. 8a
- Includes initialization routines.
- TODO: Make multi-core

DevNotes
- FIXME major refactor of this function required.
- FIXME this function actually occur during the parent CSM, therefore not all pull model #674
"""
function collectDwnInitMsgFromParent_StateMachine(csmc::CliqStateMachineContainer)
  #
  # TODO consider exit early for root clique rather than avoiding this function
  infocsm(csmc, "8a, needs down message -- attempt down init")
  setCliqDrawColor(csmc.cliq, "gold")

  # initialize clique in downward direction
  # not if parent also needs downward init message
  # prnt = getParent(csmc.tree, csmc.cliq)[1]
  opt = getSolverParams(csmc.dfg)
  @assert !haskey(opt.devParams,:dontUseParentFactorsInitDown) "dbgnew is old school, 459 dwninit consolidation has removed the option for :dontUseParentFactorsInitDown"

  # take atomic lock OF PARENT ??? when waiting for downward information
  # lockUpStatus!(prnt, prnt.index, true, csmc.logger, true, "cliq$(csmc.cliq.index)") # TODO XY ????
  # infocsm(csmc, "8a, after up lock")

  # get down message from the parent
  # check if any msgs should be multiplied together for the same variable
  # get the current messages ~~stored in~~ [going to] the parent (pull model #674)
  # FIXME, post #459 calls?
  # this guy is getting any sibling up messages by calling on the parent
  prntmsgs::Dict{Int, LikelihoodMessage} = getMsgsUpInitChildren(csmc.tree, csmc.cliq, TreeBelief, skip=[csmc.cliq.index;])         
  
  # reference to default dict location
  dwinmsgs = getfetchCliqueInitMsgDown(csmc.cliq.data, from=:getMsgDwnThisInit) |> deepcopy  #JT 459 products = getMsgDwnThisInit(prnt)
  infocsm(csmc, "getMsgInitDwnParent -- msg ids::Int=$(collect(keys(prntmsgs)))")
  
  # stack all parent incoming upward messages into dict of vector msgs
  prntBelDictVec::Dict{Symbol, Vector{TreeBelief}} = convertLikelihoodToVector(prntmsgs, logger=csmc.logger)
  ## TODO use parent factors too
  # intersect with the asking clique's separator variables
  # this function populates `dwinmsgs` with the appropriate products described in `prntBelDictVec`
  # FIXME, should not be using full .dfg ???
  condenseDownMsgsProductPrntFactors!(csmc.dfg, dwinmsgs, prntBelDictVec, prnt, csmc.cliq, csmc.logger)

  # remove msgs that have no data
  rmlist = Symbol[]
  for (prsym,belmsg) in dwinmsgs.belief
    if belmsg.inferdim < 1e-10
      # no information so remove
      push!(rmlist, prsym)
    end
  end
  infocsm(csmc, "prepCliqInitMsgsDown! -- rmlist, no inferdim, keys=$(rmlist)")
  for pr in rmlist
    delete!(dwinmsgs.belief, pr)
  end

  infocsm(csmc, "prepCliqInitMsgsDown! -- product keys=$(collect(keys(dwinmsgs.belief)))")

  # now put the newly computed message in the appropriate container
  # FIXME THIS IS A PUSH MODEL, see #674 -- must make pull model first
  # FIXME must be consolidated as part of #459
  putCliqueInitMsgDown!(getCliqueData(csmc.cliq), dwinmsgs)
  # unlock
  # unlockUpStatus!(prnt) # TODO XY ????, maybe can remove after pull model #674?
  infocsm(csmc, "8a, attemptCliqInitD., unlocked")

  # go to 7b (maybe, and part of dwnMsg #459 WIP 9)
  return slowIfChildrenNotUpSolved_StateMachine
    # # go to 8j.
    # return dwnInitSiblingWaitOrder_StateMachine
end


# currently for internal use only
# initialize variables based on best current achievable ordering
# OBVIOUSLY a lot of refactoring and consolidation needed with cliqGibbs / approxCliqMarginalUp
function initSolveSubFg!( subfg::AbstractDFG,
                          logger=ConsoleLogger() )
  #
  error("initSolveSubFg! is deprecated, do not use.  Follow post #459 CSM")
  varorder = getSubFgPriorityInitOrder(subfg, logger)
  with_logger(logger) do
    @info "initSolveSubFg! -- varorder=$varorder"
  end
  cycleInitByVarOrder!(subfg, varorder, logger=logger)
  nothing
end




# Helper function for prepCliqInitMsgsDown!
# future, be used in a cached system with parent in one location only for all siblings
# this function rebuilds a local subgraph from dfg and performs the calculations of the parent here and does not wait on the CSM to perform anything.
# 4-stroke compute may render this whole function obsolete.
function condenseDownMsgsProductPrntFactors!( fgl::AbstractDFG,
                                              products::LikelihoodMessage,
                                              msgspervar::Dict{Symbol, <:AbstractVector},
                                              prnt::TreeClique,
                                              cliq::TreeClique,
                                              logger=ConsoleLogger() )
  #
  error("condenseDownMsgsProductPrntFactors! is deprecated, follow post #459 instead")

  # determine the variables of interest
  reqMsgIds = collect(keys(msgspervar))
  # unique frontals per cliq
  prntvars = intersect(getCliqSeparatorVarIds(cliq), getCliqAllVarIds(prnt))
  lvarids = union(prntvars, reqMsgIds)
  # determine allowable factors, if any (only from parent cliq)
  awfcts = getCliqFactorIdsAll(prnt)

  # build required subgraph for parent/sibling down msgs
  lsfg = buildSubgraph(fgl, lvarids, 1; verbose=false)

  tempfcts = lsf(lsfg)
  dellist = setdiff(awfcts, tempfcts)
  for delf in dellist
    # TODO -- double check this deletefactor method is leaving the right parent sharing factor graph behind
    if exists(lsfg, delf)
      deleteFactor!(lsfg,delf)
    end
  end

  # add message priors
  addMsgFactors!(lsfg, msgspervar) # , DownwardPass

  # perform initialization/inference
  # FIXME, better consolidate with CSM ....uhhh TODO
  initSolveSubFg!(lsfg, logger)

  # extract complete downward marginal msg priors
  for id in intersect(getCliqSeparatorVarIds(cliq), lvarids)
    vari = getVariable(lsfg, id)
    products.belief[id] = TreeBelief(vari)
  end

  nothing
end
# # QUICK DBG CODE
# with_logger(logger) do
#     @info "condenseDownMsgsProductPrntFactors! -- reqMsgIds=$(reqMsgIds),"
#     @info "condenseDownMsgsProductPrntFactors! -- vars=$(lvarids),"
#     @info "condenseDownMsgsProductPrntFactors! -- allow factors $awfcts"
# end
# with_logger(logger) do
#     @info "condenseDownMsgsProductPrntFactors! -- lsfg fcts=$(lsf(lsfg)),"
#     @info "condenseDownMsgsProductPrntFactors! -- excess factors $dellist"
# end
# vars = ls(lsfg)
# len = length(vars)
# tdims = Vector{Float64}(undef, len)
# isinit = Vector{Bool}(undef, len)
# for idx in 1:len
#   tdims[idx] = getVariableSolvableDim(lsfg, vars[idx])
#   isinit[idx] = isInitialized(lsfg, vars[idx])
# end
# with_logger(logger) do
#     @info "condenseDownMsgsProductPrntFactors! -- after cycle init: vars=$vars"
#     @info "condenseDownMsgsProductPrntFactors! -- after cycle init: tdims=$tdims"
#     @info "condenseDownMsgsProductPrntFactors! -- after cycle init: isinit=$isinit"
# end
# # QUICK DBG CODE




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
  - `putMsgUpInit!(cliq, msg)`
  - `setCliqueStatus!(cliq, status)`
  - `setCliqDrawColor(cliq, "sienna")`
  - `notifyCliqDownInitStatus!(cliq, status)`

Algorithm:
- determine which downward messages influence initialization order
- initialize from singletons to most connected non-singletons
- revert back to needdownmsg if cycleInit does nothing
- can only ever return :initialized or :needdownmsg status

DevNotes
- TODO Lots of cleanup required, especially from calling function.
- TODO move directly into a CSM state function
"""
function doCliqInitDown!( subfg::AbstractDFG,
                          cliq::TreeClique,
                          initorder;
                          dbg::Bool=false,
                          logpath::String="/tmp/caesar/",
                          logger=ConsoleLogger() )
  #

  # store the cliqSubFg for later debugging
  if dbg
    DFG.saveDFG(subfg, joinpath(logpath,"logs/cliq$(cliq.index)/fg_beforedowninit"))
  end

  # cycle through vars and attempt init
  with_logger(logger) do
    @info "cliq $(cliq.index), doCliqInitDown! -- 5, cycle through vars and attempt init"
  end

  status = :needdownmsg
  if cycleInitByVarOrder!(subfg, initorder)
    status = :initialized
  end

  with_logger(logger) do
    @info "cliq $(cliq.index), doCliqInitDown! -- 6, current status: $status"
  end

  # store the cliqSubFg for later debugging
  if dbg
      DFG.saveDFG(subfg, joinpath(logpath,"logs/cliq$(cliq.index)/fg_afterdowninit"))
  end

  return status
end


function blockCliqUntilParentDownSolved(prnt::TreeClique; logger=ConsoleLogger())::Nothing
  @error("blockCliqUntilParentDownSolved is deprecated, use CSM directly")

  lbl = getLabel(prnt)

  with_logger(logger) do
    @info "blockCliqUntilParentDownSolved, prntcliq=$(prnt.index) | $lbl | going to fetch initdownchannel..."
  end
  flush(logger.stream)
  blockMsgDwnUntilStatus(prnt, :downsolved)

  return nothing
end

function blockMsgDwnUntilStatus(cliq::TreeClique, status::CliqStatus)
  @error("blockMsgDwnUntilStatus is deprecated, use CSM directly")
  while fetchDwnMsgConsolidated(cliq).status != status
    wait(getSolveCondition(cliq))
  end
  nothing
end

# """
#    $SIGNATURES
#
# Determine if this `cliq` has been fully initialized and child cliques have completed their full upward inference.
# """
# function isCliqReadyInferenceUp(fgl::FactorGraph, tree::AbstractBayesTree, cliq::TreeClique)
#   isallinit = areCliqVariablesAllInitialized(fgl, cliq)
#
#   # check that all child cliques have also completed full up inference.
#   for chl in getChildren(tree, cliq)
#     isallinit &= isUpInferenceComplete(chl)
#   end
#   return isallinit
# end

# """
#     $SIGNATURES
#
# Perform cliq initalization calculation based on current state of the tree and factor graph,
# using upward message passing logic.
#
# Notes
# - adds msg priors added to clique subgraph
# - Return either of (:initialized, :upsolved, :needdownmsg, :badinit)
# - must use factors in cliq only, ensured by using subgraph -- TODO general case.
#
# DevNotes
# - FIXME, integrate with `8f. prepInitUp_StateMachine`
# """
# function doCliqAutoInitUpPart1!(subfg::AbstractDFG,
#                                 tree::AbstractBayesTree,
#                                 cliq::TreeClique;
#                                 up_solve_if_able::Bool=true,
#                                 multiproc::Bool=true,
#                                 logger=ConsoleLogger() )
#   #
#
#   # attempt initialize if necessary
#   if !areCliqVariablesAllInitialized(subfg, cliq)
#     # structure for all up message densities computed during this initialization procedure.
#     varorder = getCliqVarInitOrderUp(tree, cliq)
#     # do physical inits, ignore cycle return value
#     cycleInitByVarOrder!(subfg, varorder, logger=logger)
#   end
#
#   return nothing
# end


"""
    $SIGNATURES

Set the marginalized status of a clique.
"""
function setCliqAsMarginalized!(cliq::TreeClique, status::Bool)
  @warn "Busy deprecating setCliqAsMarginalized!, use CSM features directly"
  if status
    setCliqueStatus!(cliq, :marginalized)
  else
    if getCliqueData(cliq).initialized == :marginalized
      @info "Reverting clique $(cliq.index) to assumed :downsolved status"
      setCliqueStatus!(cliq, :downsolved)
    else
      error("Unknown clique de-marginalization requist for clique $(cliq.index), current status: $(cliq.initialized)")
    end
  end
end

"""
    $SIGNATURES

Run through entire tree and set cliques as marginalized if all clique variables are marginalized.

Notes:
- TODO can be made fully parallel, consider converting for use with `@threads` `for`.
"""
function updateTreeCliquesAsMarginalizedFromVars!(fgl::AbstractDFG, tree::AbstractBayesTree)::Nothing
  @warn "Busy deprecating updateTreeCliquesAsMarginalizedFromVars!, use CSM features directly"
  for (clid, cliq) in getCliques(tree)
    if isCliqMarginalizedFromVars(fgl, cliq)
      setCliqAsMarginalized!(cliq, true)
    end
  end
  nothing
end


# NOTE I only saw this function after I replace all the functions with _dbgCSMSaveSubFG
# I consolidated the 2 and think this one can therefore be deprecated. Unless you use it 
# separate from a CSMC
"""
    $SIGNATURES

Internal helper function to save a dfg object to LogPath during clique state machine operations.

Notes
- will only save dfg object if `opts.dbg=true`

Related

saveDFG, loadDFG!, loadDFG
"""
function _dbgSaveDFG(dfg::AbstractDFG,
                    filename::AbstractString="fg_temp",
                    opts::AbstractParams=getSolverParams(dfg)  )::String
  #
  Base.depwarn("`_dbgSaveDFG` is deprecated use `_dbgCSMSaveSubFG`", _dbgCSMSaveSubFG)
  folder::String=joinpath(opts.logpath,"logs")
  if opts.dbg
    if !ispath(folder)
      mkpath(folder)
    end
    DFG.saveDFG(dfg, joinpath(folder, "$filename"))
    drawGraph(dfg, show=false, filepath=joinpath(folder, "$filename.pdf"))
  end
  folder*filename
end

##==============================================================================
## Delete at end v0.15.x
##==============================================================================


@deprecate putCliqueMsgDown!(cdata::BayesTreeNodeData, msg::LikelihoodMessage; from::Symbol=:nothing) putMsgDwnThis!(cdata, msg, from=from)

# """
#     $SIGNATURES

# Initialization downward message passing is different from regular inference since
# it is possible that none of the child cliq variables have been initialized.

# Notes
# - init upward msgs are individually stored in child cliques (pull model good).
# - fresh product of overlapping beliefs are calculated on each function call.
# - Assumed that `prnt` of siblings

# Dev Notes
# - This should be the initialization cycle of parent, build up bit by bit...
# """
# function prepCliqInitMsgsDown!(fgl::AbstractDFG,
#                                tree::AbstractBayesTree,
#                                prnt::TreeClique,
#                                cliq::TreeClique;
#                                logger=ConsoleLogger() )
#   #
#   # tt = split(string(now()), 'T')[end]
#   # with_logger(logger) do
#   #   @info "$(tt) prnt $(prnt.index), prepCliqInitMsgsDown! -- with cliq $(cliq.index)"
#   # end
  
#   # FIXME use only LikelihoodMessage
#   # check if any msgs should be multiplied together for the same variable
#   # msgspervar = LikelihoodMessage() # or maybe Dict{Int, LikelihoodMessage}()
#   msgspervar = getMsgInitDwnParent(tree, cliq, logger=logger)
#   # reference to default dict location
#   #JT 459 products = getMsgDwnThisInit(prnt)
#   products = getfetchCliqueInitMsgDown(prnt.data, from=:getMsgDwnThisInit) |> deepcopy
  
#   ## TODO use parent factors too
#   # intersect with the asking clique's separator variables
#   condenseDownMsgsProductPrntFactors!(fgl, products, msgspervar, prnt, cliq, logger)
  
#   # with_logger(logger) do
#   #   @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- vars fw/ down msgs=$(collect(keys(msgspervar)))"
#   # end
#   # flush(logger.stream)

#   # remove msgs that have no data
#   rmlist = Symbol[]
#   for (prsym,belmsg) in products.belief
#     if belmsg.inferdim < 1e-10
#       # no information so remove
#       push!(rmlist, prsym)
#     end
#   end
#   with_logger(logger) do
#     @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- rmlist, no inferdim, keys=$(rmlist)"
#   end
#   for pr in rmlist
#     delete!(products.belief, pr)
#   end

#   with_logger(logger) do
#     @info "cliq $(prnt.index), prepCliqInitMsgsDown! -- product keys=$(collect(keys(products.belief)))"
#   end

#   # now put the newly computed message in the appropriate container
#   # FIXME THIS IS A PUSH MODEL, see #674 -- make pull model
#   putCliqueInitMsgDown!(getCliqueData(prnt), products)

#   return products
# end


# used during nonparametric CK preparation, when information from multiple siblings must be shared together
# const IntermediateSiblingMessagesTB{T} = Vector{TreeBelief{T}}
# const IntermediateMultiSiblingMessagesTB{T} = Dict{Symbol, Vector{TreeBelief{T}}}


# FIXME, better standardize intermediate types
# can be replaced by Vector{TreeBelief}
# const IntermediateSiblingMessages = Vector{Tuple{BallTreeDensity,Float64}}
# const IntermediateMultiSiblingMessages = Dict{Symbol, IntermediateSiblingMessages}


# Helper function for prepCliqInitMsgsDown!
# populate products with products of upward messages
function condenseDownMsgsProductOnly!(fgl::AbstractDFG,
                                      products::LikelihoodMessage,
                                      msgspervar::Dict{Symbol, <:AbstractVector}  )
  #
  error("condenseDownMsgsProductOnly!(::AbstractDFG,::LikelihoodMessage, ::IntermediateMultiSiblingMessages) is obsolete")
  # multiply multiple messages together
  for (msgsym, msgsBo) in msgspervar
    # check if this particular down message requires msgsym
    if exists(fgl, msgsym) # DFG.hasVariable(fgl, msgsym)
      if length(msgspervar[msgsym]) > 1
        msgs = getindex.(msgsBo, 1)
        haspars = 0.0
        for mb in msgsBo, val in mb[2]
          haspars += val
        end
        products[msgsym] = (manifoldProduct(msgs, getManifolds(fgl, msgsym)), haspars)
      else
        # transfer if only have a single belief
        products[msgsym] = (msgsBo[1][1], msgsBo[1][2])
      end
    else
      # not required, therefore remove from message to avoid confusion
      if haskey(products, msgsym)
        delete!(products, msgsym)
      end
    end
  end
  nothing
end

@deprecate getfetchCliqueMsgDown(cdata::BayesTreeNodeData; from::Symbol=:nothing) getfetchCliqueInitMsgDown(cdata, from=from)

export getIdx

"""
    $SIGNATURES

Return interger index of desired variable element.

Example
-------
```julia
pp = RoME.Point2()
getIdx(pp, :posY) # = 2
```

Internal Notes
--------------
- uses number i < 100 for index number, and
- uses +100 offsets to track the minibatch number of the requested dimension
"""
function getIdx(pp::Tuple,
                sym::Symbol,
                i::Int=0)
  #
  error("getIdx is obsolete, use DistributedFactorGraphs objects/methods instead.")
  i-=100
  for p in pp
    i,j = getIdx(p, sym, i)
    if i > 0
      return i, j
    end
  end
  return i,-1
end
getIdx(pp::Symbol, sym::Symbol, i::Int=0) = pp==sym ? (abs(i)%100+1, div(abs(i)-100,100)) : (i-1, div(abs(i)-100,100))
function getIdx(pp::InferenceVariable, sym::Symbol, i::Int=0)
  return getIdx(pp.dimtype, sym)
end


# """
#     $SIGNATURES

# Get the latest down message from the parent node (without calculating anything).

# Notes
# - Different from down initialization messages that do calculate new values -- see `prepCliqInitMsgsDown!`.
# - Basically converts function `getDwnMsgs` from `Dict{Symbol,BallTreeDensity}` to `Dict{Symbol,Vector{BallTreeDensity}}`.

# Related

# getMsgInitDwnParent
# """
# function getMsgDwnParent(treel::AbstractBayesTree, cliq::TreeClique)
#   downmsgs = IntermediateMultiSiblingMessages()
#   prnts = getParent(treel, cliq)
#   if 0 < length(prnts)
#     prnt = prnts[1]
#     prntmsgs = getDwnMsgs(prnt)
#     for (key, bel) in prntmsgs
#       if !haskey(downmsgs, key)
#         downmsgs[key] = IntermediateSiblingMessages()
#       end
#       # TODO insert true inferred dim
#       push!(downmsgs[key], bel)
#     end
#   end

#   return downmsgs
# end


# function getMsgDwnThisInit(cdat::BayesTreeNodeData) 
#   iifdepwarn("#459 replace with getfetchCliqueMsgDown", :getMsgDwnThisInit)
#   return cdat.downInitMsg
# end

@deprecate fetchDataElement(dfg::AbstractDFG, varsym::Symbol, lbl::Symbol) fetchDataJSON(dfg, varsym, lbl)


function addMsgFactors!(subfg::AbstractDFG,
                        msgs::Dict{Symbol, <:AbstractVector} )  #Dict{Symbol, Vector{Tuple{BallTreeDensity, Float64}}} )
  # msgs::
  # add messages as priors to this sub factor graph
  @warn "Tuple{KDE,Floa64} specific version of addMsgFactors! is deprecated, use LikelihoodMessage version instead."
  msgfcts = DFGFactor[]
  svars = DFG.listVariables(subfg)
  for (msym, dms) in msgs
    for treebelief in dms
      if msym in svars
        # TODO should be on manifold prior, not just generic euclidean prior -- okay since variable on manifold, but not for long term
        fc = addFactor!(subfg, [msym], MsgPrior(manikde!(treebelief), treebelief.inferdim), graphinit=false)
        push!(msgfcts, fc)
      end
    end
  end
  return msgfcts
end



