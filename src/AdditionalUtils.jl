
"""
    $SIGNATURES

Construct and draw the subgraph associated with clique `frontalSym::Symbol`.

Notes
- See `writeGraphPdf` for details on keyword options.

Related

writeGraphPdf, buildCliqSubgraphUp, buildSubgraphFromLabels, spyCliqMat, drawTree
"""
function drawCliqSubgraphUp(fgl::FactorGraph,
                            treel::BayesTree,
                            frontalSym::Symbol;
                            show::Bool=true,
                            filepath::String="/tmp/cliq_sfg.pdf",
                            engine::AS1="sfdp",
                            viewerapp::AS2="evince"  )::FactorGraph where {AS1 <: AbstractString, AS2 <: AbstractString}
  #
  sfg = buildCliqSubgraphUp(fgl, treel, frontalSym)
  writeGraphPdf(sfg, show=show, viewerapp=viewerapp, engine=engine, filepath=filepath)
  return sfg
end
