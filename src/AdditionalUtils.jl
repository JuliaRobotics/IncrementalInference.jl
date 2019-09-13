
"""
    $SIGNATURES

Construct (new) subgraph and draw the subgraph associated with clique `frontalSym::Symbol`.

Notes
- See `drawGraphCliq`/`writeGraphPdf` for details on keyword options.

Related

drawGraphCliq, spyCliqMat, drawTree, buildCliqSubgraphUp, buildSubgraphFromLabels
"""
function drawCliqSubgraphUpMocking(fgl::G,
                                   treel::BayesTree,
                                   frontalSym::Symbol;
                                   show::Bool=true,
                                   filepath::String="/tmp/cliq_sfg.pdf",
                                   engine::AS1="sfdp",
                                   viewerapp::AS2="evince"  )::Nothing where {G <: AbstractDFG, AS1 <: AbstractString, AS2 <: AbstractString}
  #
  sfg = buildCliqSubgraphUp(fgl, treel, frontalSym)
  writeGraphPdf(sfg, show=show, viewerapp=viewerapp, engine=engine, filepath=filepath)
  nothing
end


"""
    $SIGNATURES

Draw the factor graph `fgl` as pdf and show.

Related

drawGraphCliq, drawTree
"""
function drawGraph(fgl::G;
                   viewerapp::String="evince",
                   filepath::AS="/tmp/caesar/random/fg.pdf",
                   engine::AS="neato", #sfdp
                   show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
  writeGraphPdf(fgl, filepath=filepath, show=show, viewerapp=viewerapp, engine=engine)
end

"""
    $SIGNATURES

Draw the factor graph from a clique state machine history at a particular step as pdf and show.

Related

drawCliqSubgraphUpMocking, drawGraph, drawTree
"""
function drawGraphCliq(hists::Dict{Int, <: Tuple},
                       step::Int,
                       tree::BayesTree,
                       frontal::Symbol;
                       show::Bool=true  )
  #
  cid = getCliq(tree, frontal).index
  cfg = hists[cid][step][4].cliqSubFg
  writeGraphPdf(cfg, show=show)
end
