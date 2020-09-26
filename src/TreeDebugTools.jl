


"""
    $(SIGNATURES)

Calculate a fresh (single step) approximation to the variable `sym` in clique `cliq` as though during the upward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variables based on the structure of the clique, factors, and incoming messages.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `getClique(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductUp(fg::AbstractDFG,
                        tree::AbstractBayesTree,
                        cliq::Symbol,
                        sym::Symbol;
                        N::Int=100,
                        dbg::Bool=false  )
  #
  cliq = getClique(tree, cliq)
  cliqdata = getCliqueData(cliq)

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  # upmsgssym = LikelihoodMessage[]
  # for cl in childCliqs(tree, cliq)
  #   msgdict = getUpMsgs(cl)
  #   dict = Dict{Symbol, TreeBelief}()
  #   for (dsy, btd) in msgdict.belief
  #     vari = getVariable(fg, dsy)
  #     # manis = getSofttype(vari).manifolds
  #     dict[dsy] = TreeBelief(btd.val, btd.bw, btd.inferdim, getSofttype(vari))
  #   end
  #   push!( upmsgssym, LikelihoodMessage(beliefDict=dict) )
  # end
  upmsgssym = getMsgsUpChildren(tree, cliq, TreeBelief)

  # perform the actual computation
  manis = getSofttype(getVariable(fg, sym)) |> getManifolds
  pGM, potprod, fulldim = cliqGibbs( fg, cliq, sym, upmsgssym, N, dbg, manis )

  return pGM, potprod
end



"""
    $(SIGNATURES)

Calculate a fresh---single step---approximation to the variable `sym` in clique `cliq` as though during the downward message passing.  The full inference algorithm may repeatedly calculate successive apprimxations to the variable based on the structure of variables, factors, and incoming messages to this clique.
Which clique to be used is defined by frontal variable symbols (`cliq` in this case) -- see `getClique(...)` for more details.  The `sym` symbol indicates which symbol of this clique to be calculated.  **Note** that the `sym` variable must appear in the clique where `cliq` is a frontal variable.
"""
function treeProductDwn(fg::G,
                        tree::AbstractBayesTree,
                        cliq::Symbol,
                        sym::Symbol;
                        N::Int=100,
                        dbg::Bool=false  ) where G <: AbstractDFG
  #
  @warn "treeProductDwn might not be working properly at this time. (post DFG v0.6 upgrade maintenance required)"
  cliq = getClique(tree, cliq)
  cliqdata = getCliqueData(cliq)

  # get the local variable id::Int identifier
  # vertid = fg.labelDict[sym]

  # get all the incoming (upward) messages from the tree cliques
  # convert incoming messages to Int indexed format (semi-legacy format)
  cl = parentCliq(tree, cliq)
  msgdict = getDwnMsgs(cl[1])
  dict = Dict{Int, TreeBelief}()
  for (dsy, btd) in msgdict
      dict[fg.IDs[dsy]] = TreeBelief(btd.val, btd.bw, btd.inferdim, getSofttype(getVariable(fg,sym)) )
  end
  dwnmsgssym = LikelihoodMessage[LikelihoodMessage(dict);]

  # perform the actual computation
  pGM, potprod, fulldim = cliqGibbs( fg, cliq, sym, dwnmsgssym, N, dbg ) #vertid

  return pGM, potprod, sym, dwnmsgssym
end