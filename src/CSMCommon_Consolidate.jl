


"""
    $SIGNATURES

Special internal function to try return the clique data if succesfully identified in `othertree::AbstractBayesTree`,
based on contents of `seeksSimilar::BayesTreeNodeData`.

Notes
- Used to identify and skip similar cliques (i.e. recycle computations)
"""
function attemptTreeSimilarClique(othertree::AbstractBayesTree, 
                                  seeksSimilar::BayesTreeNodeData  )
  #
  # inner convenience function for returning empty clique
  function EMPTYCLIQ()
    clq = TreeClique(-1,"null")
    setLabel!(clq, "")
    setCliqueData!(clq, BayesTreeNodeData())
    return clq
  end

  # does the other clique even exist?
  seekFrontals = getCliqFrontalVarIds(seeksSimilar)
  if !hasClique(othertree, seekFrontals[1])
    return EMPTYCLIQ()
  end

  # do the cliques share the same frontals?
  otherCliq = getClique(othertree, seekFrontals[1])
  otherFrontals = getCliqFrontalVarIds(otherCliq)
  commonFrontals = intersect(seekFrontals, otherFrontals)
  if length(commonFrontals) != length(seekFrontals) || length(commonFrontals) != length(otherFrontals)
    return EMPTYCLIQ()
  end

  # do the cliques share the same separator variables?
  seekSeparator = getCliqSeparatorVarIds(seeksSimilar)
  otherSeparator = getCliqSeparatorVarIds(otherCliq)
  commonSep = intersect(seekSeparator, otherSeparator)
  if length(commonSep) != length(seekSeparator) || length(commonSep) != length(otherSeparator)
    return EMPTYCLIQ()
  end

  # do the cliques use the same factors (potentials)
  seekPotentials = getCliqFactorIds(seeksSimilar)
  otherFactors = getCliqFactorIds(otherCliq)
  commonFactors = intersect(seekPotentials, otherFactors)
  if length(commonFactors) != length(seekPotentials) || length(commonFactors) != length(otherFactors)
    return EMPTYCLIQ()
  end

  # lets assume they are the same
  return otherCliq::TreeClique
end
