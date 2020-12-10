export getVariableOrder, calcCliquesRecycled
export getCliquePotentials
export getClique, getCliques, getCliqueIds, getCliqueData
export hasClique
export setCliqueDrawColor!, getCliqueDrawColor

"""
    $SIGNATURES

Set the color of a cliq in the Bayes (Junction) tree.
"""
function setCliqueDrawColor!(cliq::TreeClique, fillcolor::String)::Nothing
  cliq.attributes["fillcolor"] = fillcolor
  cliq.attributes["style"] = "filled"
  nothing
end


function getCliqueDrawColor(cliq::TreeClique)
  haskey(cliq.attributes, "fillcolor") ? cliq.attributes["fillcolor"] : nothing
end

"""
    $SIGNATURES

Get the frontal variable IDs `::Int` for a given clique in a Bayes (Junction) tree.
"""
getCliqFrontalVarIds(cliqdata::Union{BayesTreeNodeData, PackedBayesTreeNodeData})::Vector{Symbol} = cliqdata.frontalIDs
getCliqFrontalVarIds(cliq::TreeClique)::Vector{Symbol} = getCliqFrontalVarIds(getCliqueData(cliq))

"""
    $SIGNATURES

Get the frontal variable IDs `::Int` for a given clique in a Bayes (Junction) tree.
"""
getFrontals(cliqd::Union{TreeClique,BayesTreeNodeData})::Vector{Symbol} = getCliqFrontalVarIds(cliqd)


"""
    $SIGNATURES

Create a new clique.
"""
function addClique!(bt::MetaBayesTree,
                    dfg::AbstractDFG,
                    varID::Symbol,
                    condIDs::Array{Symbol}=Symbol[]  )
  bt.btid += 1 #increment cliqueID counter

  cId = CliqueId(bt.btid)
  
  clq = TreeClique(cId)
  setLabel!(clq, "")

  #TODO addClique!(bt, clq), can't we already have the parent here
  # if isa(bt.bt,GenericIncidenceList)
  #   Graphs.add_vertex!(bt.bt, clq)
  #   bt.cliques[bt.btid] = clq
  #   clId = bt.btid
  # if isa(bt.bt, MetaDiGraph)
  #   @assert MetaGraphs.add_vertex!(bt.bt, :clique, clq) "add_vertex! failed"
  #   clId = MetaGraphs.nv(bt.bt)
  #   MetaGraphs.set_indexing_prop!(bt.bt, clId, :cliqId, bt.btid)
  # else
  #   error("Oops, something went wrong when adding a clique to the tree")
  # end
  @assert MetaGraphs.add_vertex!(bt.bt, :clique, clq) "Error trying to addClique! - add_vertex! failed"
  MetaGraphs.set_indexing_prop!(bt.bt, MetaGraphs.nv(bt.bt), :cliqId, cId)

  appendClique!(bt, cId, dfg, varID, condIDs)
  return clq
end


"""
    $(SIGNATURES)

Return the TreeClique node object that represents a clique in the Bayes
(Junction) tree, as defined by one of the frontal variables `frt<:AbstractString`.

Notes
- Frontal variables only occur once in a clique per tree, therefore is a unique identifier.

Related:

getCliq, getTreeAllFrontalSyms
"""
getClique(tree::AbstractBayesTree, cId::Int) = tree.cliques[cId]
getClique(bt::AbstractBayesTree, frt::Symbol) = getClique(bt, bt.frontals[frt])
getClique(tree::MetaBayesTree, cId::CliqueId) = MetaGraphs.get_prop(tree.bt, tree.bt[:cliqId][cId], :clique)
getClique(tree::MetaBayesTree, cIndex::Int) = MetaGraphs.get_prop(tree.bt, cIndex, :clique)

"""
    $(SIGNATURES)
"""
function deleteClique!(tree::MetaBayesTree, cId::CliqueId)
  clique = getClique(tree, cId) 
  @assert MetaGraphs.rem_vertex!(tree.bt, tree.bt[:cliqId][cId]) "rem_vertex! failed"
  foreach(frt->delete!(tree.frontals,frt), getFrontals(clique))
  return clique
end

isRoot(tree::MetaBayesTree, cliq::TreeClique) = isRoot(tree, cliq.id)
function isRoot(tree::MetaBayesTree, cliqId::CliqueId)
  length(MetaGraphs.inneighbors(tree.bt, tree.bt[:cliqId][cliqId])) == 0
end

"""
    $SIGNATURES

Return boolean on whether the frontal variable `frt::Symbol` exists somewhere in the `::BayesTree`.
"""
hasClique(bt::AbstractBayesTree, frt::Symbol) = haskey(bt.frontals, frt)

"""
    $SIGNATURES

Return depth in tree as `::Int`, with root as depth=0.

Related

getCliq
"""
function getCliqDepth(tree, cliq)::Int
  prnt = getParent(tree, cliq)
  if length(prnt) == 0
    return 0
  end
  return getCliqDepth(tree, prnt[1]) + 1
end
getCliqDepth(tree::AbstractBayesTree, sym::Symbol)::Int = getCliqDepth(tree, getClique(tree, sym))


getCliques(tree::AbstractBayesTree) = tree.cliques

function getCliques(tree::MetaBayesTree)
  d = Dict{Int,Any}()
  for (k,v) in tree.bt.vprops
    d[k] = v[:clique]
  end
  return d
end

getCliqueIds(tree::AbstractBayesTree) = keys(getCliques(tree))

function getCliqueIds(tree::MetaBayesTree)
  # MetaGraphs.vertices(tree.bt)
  keys(tree.bt.metaindex[:cliqId])
end

"""
    $SIGNATURES

Return reference to the clique data container.
"""
getCliqueData(cliq::TreeClique) = cliq.data
getCliqueData(tree::AbstractBayesTree, cId::Union{CliqueId, Int}) = getClique(tree, cId) |> getCliqueData

"""
    $SIGNATURES

Set the clique data container to a new object `data`.
"""
setCliqueData!(cliq::TreeClique, data::Union{PackedBayesTreeNodeData, BayesTreeNodeData}) = cliq.data = data
setCliqueData!(tree::AbstractBayesTree, cId::Int, data::Union{PackedBayesTreeNodeData, BayesTreeNodeData}) = setCliqueData!(getClique(tree, cId), data)

# #TODO
# addClique!(tree::AbstractBayesTree, parentCliqId::Int, cliq::TreeClique) = error("addClique!(tree::AbstractBayesTree, parentCliqId::Int, cliq::TreeClique) not implemented")
# updateClique!(tree::AbstractBayesTree, cliq::TreeClique) = error("updateClique!(tree::AbstractBayesTree, cliq::TreeClique)::Bool not implemented")
# deleteClique!(tree::AbstractBayesTree, cId::Int) =  error("deleteClique!(tree::AbstractBayesTree, cId::Int)::TreeClique not implemented")

"""
    $SIGNATURES

Generate the label for particular clique (used by graphviz for visualization).
"""
function makeCliqueLabel(dfg::G, bt::AbstractBayesTree, clqID::CliqueId)::String where G <: AbstractDFG
  clq = getClique(bt, clqID)
  flbl = ""
  clbl = ""
  for fr in getCliqueData(clq).frontalIDs
    flbl = string(flbl, DFG.getVariable(dfg,fr).label, ",")
  end
  for sepr in getCliqueData(clq).separatorIDs
    clbl = string(clbl, DFG.getVariable(dfg,sepr).label, ",")
  end
  setLabel!(clq, string(clqID,"| ", flbl, ": ", clbl))
end

"""
    $SIGNATURES

Add the separator for the newly created clique.
"""
function appendSeparatorToClique(bt::AbstractBayesTree, clqID::CliqueId, seprIDs::Array{Symbol,1})
  union!(getCliqueData(bt, clqID).separatorIDs, seprIDs)
  nothing
end

"""
    $SIGNATURES

Add a new frontal variable to clique.

DevNotes
- TODO, define what "conditionals" are CLEARLY!!
"""
function appendClique!(bt::AbstractBayesTree,
                       clqID::CliqueId,
                       dfg::AbstractDFG,
                       varID::Symbol,
                       seprIDs::Array{Symbol,1}=Symbol[] )::Nothing
  #
  clq = getClique(bt, clqID)
  var = DFG.getVariable(dfg, varID)

  # add frontal variable
  push!(getCliqueData(clq).frontalIDs, varID)

  # total dictionary of frontals for easy access
  bt.frontals[varID] = clqID

  # TODO - confirm this, append to cliq separator ??
  # @info "going for: appendSeparatorToClique on (clqID, seprIDs)=($clqID, $seprIDs)"
  appendSeparatorToClique(bt, clqID, seprIDs)
  makeCliqueLabel(dfg, bt, clqID)
  return nothing
end

"""
    $SIGNATURES

Instantiate a new child clique in the tree.
"""
function newChildClique!(bt::AbstractBayesTree, dfg::AbstractDFG, CpID::CliqueId, varID::Symbol, Sepj::Array{Symbol,1})
  #TODO this is too manual, replace with just addClique that takes as argument the parent
        # addClique!(bt, dfg, CpID, varID, Sepj)
  # physically create the new clique
  chclq = addClique!(bt, dfg, varID, Sepj)
  parent = getClique(bt, CpID)

  # if isa(bt.bt,GenericIncidenceList)
  #   # Staying with Graphs.jl for tree in first stage
  #   edge = Graphs.make_edge(bt.bt, parent, chclq)
  #   Graphs.add_edge!(bt.bt, edge)
  if isa(bt.bt, MetaDiGraph)
    # TODO EDGE properties here
    @assert MetaGraphs.add_edge!(bt.bt, bt.bt[:cliqId][CpID], bt.bt[:cliqId][chclq.id]) "Add edge failed"
  else
    error("Oops, something went wrong when adding a new child clique to the tree")
  end


  return chclq
end

"""
    $SIGNATURES

Post order tree traversal and build potential functions
"""
function findCliqueFromFrontal(bt::AbstractBayesTree, frtlID::Int)
  for cliqPair in getCliques(bt, cliques)
    id = cliqPair[1]
    cliq = cliqPair[2]
    for frtl in getFrontals(cliq)
      if frtl == frtlID
        return cliq
      end
    end
  end
  error("Clique with desired frontal ID not found")
end

"""
    $SIGNATURES

Find parent clique Cp that containts the first eliminated variable of Sj as frontal.
"""
function identifyFirstEliminatedSeparator(dfg::AbstractDFG,
                                          elimorder::Vector{Symbol},
                                          firvert::DFGVariable,
                                          Sj=getSolverData(firvert).separator)::DFGVariable
  #
  firstelim = 99999999999
  for s in Sj
    temp = something(findfirst(isequal(s), elimorder), 0) # findfirst(p, s)
    if (temp < firstelim)
      firstelim = temp
    end
  end
  DFG.getVariable(dfg, elimorder[firstelim])
end

"""
    $SIGNATURES

Eliminate a variable and add to tree cliques accordingly.

Dev Notes
- `p` should be elimination order.
- `var` is next variable to be added to the tree.
- TODO, make sure this works for disjoint graphs
  - Check laplacian, check eigen == 1 is disjoint sets

References
Kaess et al.: Bayes Tree, WAFR, 2010, [Alg. 2]
Kaess et al.: iSAM2, IJRR, 2011, [Alg. 3]
Fourie, D.: mmisam, PhD thesis, 2017. [Chpt. 5]
"""
function newPotential(tree::AbstractBayesTree, dfg::G, var::Symbol, elimorder::Array{Symbol,1}) where G <: AbstractDFG
    firvert = DFG.getVariable(dfg,var)
    # no parent
    if (length(getSolverData(firvert).separator) == 0)
      # if (length(getCliques(tree)) == 0)
        # create new root
        addClique!(tree, dfg, var)
      # else
      #   # add to root
      #   @warn "root append clique is happening"
      #   appendClique!(tree, 1, dfg, var)
      # end
    else
      # find parent clique Cp that containts the first eliminated variable of Sj as frontal
      Sj = getSolverData(firvert).separator
      felbl = identifyFirstEliminatedSeparator(dfg, elimorder, firvert, Sj).label
      # get clique id of first eliminated frontal
      CpID = tree.frontals[felbl]
      # look to add this conditional to the tree
      cliq = getClique(tree, CpID)
      # clique of the first eliminated frontal
      unFC = union(getCliqFrontalVarIds(cliq), getCliqSeparatorVarIds(cliq))
      # if the separator of this new variable is identical to the (entire) clique of the firstly eliminated frontal.
      if (sort(unFC) == sort(Sj))
        # just add new variable as frontal to this clique
        # insert conditional (p(var|sepr)) into clique CpID -- i.e. just adding a frontal
        # @info "adding new frontal $var to existing clique $CpID"
        appendClique!(tree, CpID, dfg, var)
      else
        # a new child clique is required here (this becomes parent)
        # @info "adding new child clique with parent separator."
        newChildClique!(tree, dfg, CpID, var, Sj)
      end
    end
end

"""
    $SIGNATURES

Build the whole tree in batch format.
"""
function buildTree!(tree::AbstractBayesTree, dfg::AbstractDFG, elimorder::Array{Symbol,1})
  revorder = reverse(elimorder,dims=1) # flipdim(p, 1), fixing #499
  # prevVar = 0
  for var in revorder
    @debug "Adding $var to tree..."
    newPotential(tree, dfg, var, elimorder)
    prevVar = var
  end

  return tree
end

"""
    $SIGNATURES

Open view to see the graphviz exported Bayes tree, assuming default location and
viewer app.  See keyword arguments for more details.
"""
function showTree(; filepath::String="/tmp/caesar/bt.dot",
                    viewerapp::String="xdot"  )
  #
  try
    @async run(`$(viewerapp) $(filepath)`)
  catch ex
    @warn "not able to show via $(viewerapp) $(filepath)"
    @show ex
    @show stacktrace()
  end
end


# A replacement for _to_dot that saves only plotting attributes
function savedot_attributes(io::IO, g::MetaDiGraph)
    write(io, "digraph G {\n")
    for p in props(g)
        write(io, "$(p[1])=$(p[2]);\n")
    end

    for v in MetaGraphs.vertices(g)
        write(io, "$v")
        if length(props(g, v)) > 0
            write(io, " [ ")
        end
        for p in props(g, v)
            # key = p[1]
            # write(io, "$key=\"$(p[2])\",")
            for (k,v) in p[2]
              write(io, "\"$k\"=\"$v\",")
            end
        end
        if length(props(g, v)) > 0
            write(io, "];")
        end
        write(io, "\n")
    end

    for e in MetaGraphs.edges(g)
        write(io, "$(MetaGraphs.src(e)) -> $(MetaGraphs.dst(e)) [ ")
        if MetaGraphs.has_prop(g, e, :downMsg) && MetaGraphs.has_prop(g, e, :upMsg)
          if isready(MetaGraphs.get_prop(g, e, :downMsg))
            write(io, "color=red")
          elseif isready(MetaGraphs.get_prop(g, e, :upMsg))
            write(io, "color=orange")
          else
            write(io, "color=black")
          end
        end
        write(io, "]\n")
    end
    write(io, "}\n")
end

function _to_dot(mdigraph::MetaDiGraph)
  g = deepcopy(mdigraph)
  for (i,val) in g.vprops
    push!(g.vprops[i],:attributes=>val[:clique].attributes)
    delete!(g.vprops[i],:clique)
    delete!(g.vprops[i],:cliqId)
  end
  m = PipeBuffer()
  savedot_attributes(m, g)
  data = take!(m)
  close(m)
  return String(data)
end


"""
    $SIGNATURES

Draw the Bayes (Junction) tree by means of `.dot` files and `pdf` reader.

Notes
- Uses system install of graphviz.org.
- Can also use Linux tool `xdot`.
- `xlabels` is optional `cliqid=>xlabel`.
"""
function drawTree(treel::AbstractBayesTree;
                  show::Bool=false,                  # must remain false for stability and automated use in solver
                  suffix::AbstractString="_"*(split(string(uuid1()),'-')[1]),
                  filepath::String="/tmp/caesar/random/bt$(suffix).dot",
                  xlabels::Dict{Int,String}=Dict{Int,String}(),
                  dpi::Int=200,
                  viewerapp::String="xdot",
                  imgs::Bool=false )
  #
  fext = filepath[(end-2):end] #split(filepath, '.')[end]
  fpwoext = filepath[1:(end-4)]# split(filepath, '.')[end-1]
  path = mkpath(dirname(fpwoext))
  
  # modify a deepcopy
  btc = deepcopy(treel)
  for (cid, cliq) in getCliques(btc)
    if imgs
      firstlabel = split(getLabel(cliq),',')[1]
      spyCliqMat(cliq, suppressprint=true) |> exportimg(path*"/$firstlabel$suffix.png")
      cliq.attributes["image"] = path*"/$firstlabel$suffix.png"
      setLabel!(cliq, "")
    end
    delete!(cliq.attributes, "data")
  end
  
  # add any xlabel info
  for (cliqid, xlabel) in xlabels
    btc.bt.vertices[cliqid].attributes["xlabel"] = xlabel
  end

  fid = IOStream("")
  try
    fid = open("$(fpwoext).dot","w+")
    write(fid,_to_dot(btc.bt))
    close(fid)
    if string(fext) == "png"
      run(`dot $(fpwoext).dot -T $(fext) -Gdpi=$dpi -o $(filepath)`)
    else
      run(`dot $(fpwoext).dot -T $(fext) -o $(filepath)`)
    end
  catch ex
    @warn ex
    @show stacktrace()
  finally
    close(fid)
  end

  show ? showTree(viewerapp=viewerapp, filepath=filepath) : nothing
end

"""
    $SIGNATURES

If opt.drawtree then start an async task to draw tree in a loop according to opt.drawtreerate.

Notes
- wont draw if opt.drawtree=false, just skips back to caller.
- Currently @async
- use `opt.showtree::Bool`
- Does not work too well when opt.async during solveTree! call, but user can use this function separately.

DevNotes
- TODO, use Threads.@spawn instead.

Related

drawTree, drawGraph
"""
function drawTreeAsyncLoop( tree::AbstractBayesTree,
                            opt::SolverParams;
                            filepath=joinLogPath(opt,"bt.pdf"),
                            dotreedraw = Int[1;]  )
  #
  # single drawtreerate
  treetask = if opt.drawtree
    @async begin
      xlabels = Dict{Int,String}()
      while dotreedraw[1] == 1 && 0 < opt.drawtreerate
        # actually draw the tree
        drawTree(tree,show=false,filepath=filepath)
        sleep(1/opt.drawtreerate)
      end
      drawTree(tree,show=opt.showtree,filepath=filepath)
    end
  end
  return treetask, dotreedraw
end

"""
    $SIGNATURES

Draw the Bayes (junction) tree with LaTeX labels by means of `.dot` and `.tex`
files.

Notes
- Uses system install of graphviz.org.
- Uses external python `dot2tex` tool (`pip install dot2tex`).

Related:

drawTree
"""
function generateTexTree( treel::AbstractBayesTree;
                          filepath::String="/tmp/caesar/bayes/bt")
  #
  btc = deepcopy(treel)
  for (cid, cliq) in getCliques(btc)
    label = getLabel(cliq)

    # Get frontals and separator, and split into elements.
    frt, sep = split(label,':')
    efrt = split(frt,',')
    esep = split(sep,',')

    # Transform frontals into latex.
    newfrontals = ""
    for l in efrt
      # Split into symbol and subindex (letter and number).
      parts = split(l, r"[^a-z0-9]+|(?<=[a-z])(?=[0-9])|(?<=[0-9])(?=[a-z])")
      if size(parts)[1] == 2
        newfrontals = string(newfrontals, "\\bm{", parts[1], "}_{", parts[2], "}, ")
      elseif size(parts)[1] == 3
        newfrontals = string(newfrontals, "\\bm{", parts[2], "}_{", parts[3], "}, ")
      end
    end
    # Remove the trailing comma.
    newfrontals = newfrontals[1:end-2]

    # Transform separator into latex.
    newseparator = ""
    if length(sep) > 1
      for l in esep
        # Split into symbol and subindex.
        parts = split(l, r"[^a-z0-9]+|(?<=[a-z])(?=[0-9])|(?<=[0-9])(?=[a-z])")
        if size(parts)[1] == 2
          newseparator = string(newseparator, "\\bm{", parts[1], "}_{", parts[2], "}, ")
        elseif size(parts)[1] == 3
          newseparator = string(newseparator, "\\bm{", parts[2], "}_{", parts[3], "}, ")
        end
      end
    end
    # Remove the trailing comma.
    newseparator = newseparator[1:end-2]
    # Create full label and replace the old one.
    newlabel = string(newfrontals, ":", newseparator)
    setLabel!(cliq, newlabel)
  end

  # Use new labels to produce `.dot` and `.tex` files.
  fid = IOStream("")
  try
    mkpath(joinpath((split(filepath,'/')[1:(end-1)])...) )
    fid = open("$(filepath).dot","w+")
    write(fid, _to_dot(btc.bt))
    close(fid)
    # All in one command.
    run(`dot2tex -tmath --preproc $(filepath).dot -o $(filepath)proc.dot`)
    run(`dot2tex $(filepath)proc.dot -o $(filepath).tex`)
  catch ex
    @warn ex
    @show stacktrace()
  finally
    close(fid)
  end

  return btc
end

"""
    $SIGNATURES

Build Bayes/Junction/Elimination tree from a given variable ordering.
"""
function buildTreeFromOrdering!(dfg::InMemoryDFGTypes,
                                p::Vector{Symbol};
                                drawbayesnet::Bool=false,
                                maxparallel::Union{Nothing, Int}=nothing,
                                solvable::Int=1 )
  #

  t0 =time_ns()
  println()
  fge = deepcopy(dfg)
  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    getSolverParams(fge).maxincidence = maxparallel
  end
  @info "Building Bayes net..."
  buildBayesNet!(fge, p, solvable=solvable)

  @info "Staring the Bayes tree construction from Bayes net"
  tree = emptyBayesTree()
  tree.variableOrder = p
  buildTree!(tree, fge, p)

  if drawbayesnet
    println("Bayes Net")
    sleep(0.1)
    fid = open("bn.dot","w+")
    write(fid,_to_dot(fge.bn))
    close(fid)
  end

  @debug "Find potential functions for each clique"
  for cliqIds in getCliqueIds(tree)
    if isRoot(tree, cliqIds)
      cliq = getClique(tree, cliqIds) # start at the root
      buildCliquePotentials(dfg, tree, cliq, solvable=solvable); # fg does not have the marginals as fge does
    end
  end
  tree.buildTime = (time_ns()-t0)/1e9
  return tree
end


"""
    $SIGNATURES

Build Bayes/Junction/Elimination tree from a given variable ordering.
"""
function buildTreeFromOrdering!(dfg::DFG.AbstractDFG,
                                p::Vector{Symbol};
                                drawbayesnet::Bool=false,
                                maxparallel::Union{Nothing, Int}=nothing  )
  #
  @debug "Copying to a local DFG"
  fge = InMemDFGType(solverParams=getSolverParams(dfg))
    #TODO JT - I think an optional solvable filter is needed in buildTreeFromOrdering!
  DFG.deepcopyGraph!(fge, dfg)
  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    getSolverParams(fge).maxincidence = maxparallel
  end

  println("Building Bayes net from cloud...")
  buildBayesNet!(fge, p)

  tree = emptyBayesTree()
  tree.variableOrder = p
  buildTree!(tree, fge, p)

  if drawbayesnet
    println("Bayes Net")
    sleep(0.1)
    fid = open("bn.dot","w+")
    write(fid,_to_dot(fge.bn))
    close(fid)
  end

  println("Find potential functions for each clique")
  cliq = getClique(tree, 1) # start at the root
  buildCliquePotentials(dfg, tree, cliq); # fg does not have the marginals as fge does

  return tree
end


"""
    $SIGNATURES

Build Bayes/Junction/Elimination tree.

Notes
- Default to free qr factorization for variable elimination order.
"""
function prepBatchTree!(dfg::AbstractDFG;
                        variableOrder::Union{Nothing, Vector{Symbol}}=nothing,
                        variableConstraints::Vector{Symbol}=Symbol[],
                        ordering::Symbol= 0==length(variableConstraints) ? :qr : :ccolamd,
                        drawpdf::Bool=false,
                        show::Bool=false,
                        filepath::String="/tmp/caesar/bt.pdf",
                        viewerapp::String="evince",
                        imgs::Bool=false,
#                         drawbayesnet::Bool=false,
                        maxparallel::Union{Nothing, Int}=nothing )
  #
  p = variableOrder !== nothing ? variableOrder : getEliminationOrder(dfg, ordering=ordering, constraints=variableConstraints)

  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    getSolverParams(dfg).maxincidence = maxparallel
  end

  # for debuggin , its useful to have the variable ordering
  if drawpdf
    ispath(getLogPath(dfg)) ? nothing : Base.mkpath(getLogPath(dfg))
    open(joinLogPath(dfg,"variableOrder.txt"), "a") do io
      writedlm(io, string.(reshape(p,1,:)), ',')
    end
  end

  tree = buildTreeFromOrdering!(dfg, p, drawbayesnet=false) # drawbayesnet

  # GraphViz.Graph(_to_dot(tree.bt))
  # Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
  #Michael reference 3sig -- x2l1x4x3    x1|x2
  @info "Bayes Tree Complete"
  if drawpdf
    drawTree(tree, show=show, filepath=filepath, viewerapp=viewerapp, imgs=imgs)
  end

  return tree
end

"""
    $SIGNATURES

Partial reset of basic data fields in `::VariableNodeData` of `::FunctionNode` structures.
"""
function resetData!(vdata::VariableNodeData)::Nothing
  vdata.eliminated = false
  vdata.BayesNetOutVertIDs = Symbol[]
  # vdata.BayesNetVertID = :_null # TODO dont use nothing, see DFG issue #16
  vdata.separator = Symbol[]
  nothing
end

function resetData!(vdata::FunctionNodeData)::Nothing
  vdata.eliminated = false
  vdata.potentialused = false
  nothing
end

"""
    $SIGNATURES

Wipe data from `dfg` object so that a completely fresh Bayes/Junction/Elimination tree
can be constructed.
"""
function resetFactorGraphNewTree!(dfg::AbstractDFG)
  for v in DFG.getVariables(dfg)
    resetData!(getSolverData(v))
  end
  for f in DFG.getFactors(dfg)
    resetData!(getSolverData(f))
  end
  nothing
end

"""
    $SIGNATURES

Reset factor graph and build a new tree from the provided variable ordering `p`.

Related

resetBuildTree!
"""
function resetBuildTreeFromOrder!(fgl::AbstractDFG, p::Vector{Symbol})
  resetFactorGraphNewTree!(fgl)
  return buildTreeFromOrdering!(fgl, p)
end

"""
    $(SIGNATURES)

Build a completely new Bayes (Junction) tree, after first wiping clean all
temporary state in fg from a possibly pre-existing tree.

Related:

buildTreeFromOrdering!, resetBuildTreeFromOrder!
"""
function resetBuildTree!( dfg::AbstractDFG;
                          ordering::Symbol=:qr,
                          drawpdf::Bool=false,
                          show::Bool=false,
                          filepath::String="/tmp/caesar/bt.pdf",
                          viewerapp::String="evince",
                          imgs::Bool=false,
                          maxparallel::Union{Nothing, Int}=nothing,
                          ensureSolvable::Bool=true,
                          variableOrder::Union{Nothing, Vector{Symbol}}=nothing,
                          variableConstraints::Vector{Symbol}=Symbol[]  )
  #
  if ensureSolvable
    ensureSolvable!(dfg)
  end
  # depcrecation
  if maxparallel !== nothing
    @warn "maxparallel keyword is deprecated, use getSolverParams(fg).maxincidence instead."
    getSolverParams(dfg).maxincidence = maxparallel
  end

  resetFactorGraphNewTree!(dfg);
  return prepBatchTree!(dfg, variableOrder=variableOrder, ordering=ordering, drawpdf=drawpdf, show=show, filepath=filepath, viewerapp=viewerapp, imgs=imgs, variableConstraints=variableConstraints);
end


"""
    $(SIGNATURES)
Experimental create and initialize tree message channels
"""
function initTreeMessageChannels!(tree::MetaBayesTree)
  for e = MetaGraphs.edges(tree.bt)
    set_props!(tree.bt, e, Dict{Symbol,Any}(:upMsg=>Channel{LikelihoodMessage}(0),:downMsg=>Channel{LikelihoodMessage}(0)))
    # push!(tree.messageChannels, e=>(upMsg=Channel{LikelihoodMessage}(0),downMsg=Channel{LikelihoodMessage}(0)))
  end
  return nothing
end

function appendUseFcts!(usefcts,
                        lblid::Symbol,
                        fct::DFGFactor )
                        # fid::Symbol )
  #
  union!(usefcts, Symbol(fct.label))
  nothing
end


"""
    $SIGNATURES

Get all factors connected to frontal variables

Dev Notes
- Why not just do this, `ffcs = union( map(x->ls(fgl, x), frtl) )`?
"""
function getCliqFactorsFromFrontals(fgl::G,
                                    cliq::TreeClique,
                                    varlist::Vector{Symbol};
                                    inseparator::Bool=true,
                                    unused::Bool=true,
                                    solvable::Int=1  ) where {G <: AbstractDFG}
  #
  frtls = getCliqueData(cliq).frontalIDs
  seprs = getCliqueData(cliq).separatorIDs
  allids = [frtls;seprs]

  usefcts = Symbol[]
  for frsym in frtls
    # usefcts = Int[]
    for fctid in ls(fgl, frsym)
        fct = getFactor(fgl, fctid)
        if !unused || !getSolverData(fct).potentialused
            loutn = ls(fgl, fctid, solvable=solvable)
            # deal with unary factors
            if length(loutn)==1
                union!(usefcts, Symbol[Symbol(fct.label);])
                # appendUseFcts!(usefcts, loutn, fct) # , frsym)
                getSolverData(fct).potentialused = true
            end
            # deal with n-ary factors
            for sep in loutn
                if sep == frsym
                    continue # skip the frsym itself
                end
                insep = sep in allids
                if !inseparator || insep
                    union!(usefcts, Symbol[Symbol(fct.label);])
                    getSolverData(fct).potentialused = true
                    if !insep
                        @debug "cliq=$(cliq.id) adding factor that is not in separator, $sep"
                    end
                end
            end
        end
    end
  end

  return usefcts
end

"""
    $SIGNATURES

Return `::Bool` on whether factor is a partial constraint.
"""
isPartial(fcf::T) where {T <: FunctorInferenceType} = :partial in fieldnames(T)
function isPartial(fct::DFGFactor)  #fct::TreeClique
  fcf = getSolverData(fct).fnc.usrfnc!
  isPartial(fcf)
end

"""
    $SIGNATURES

Determine and set the potentials for a particular `cliq` in the Bayes (Junction) tree.
"""
function setCliqPotentials!(dfg::G,
                            bt::AbstractBayesTree,
                            cliq::TreeClique;
                            solvable::Int=1  ) where G <: AbstractDFG
  #
  varlist = getCliqVarIdsAll(cliq)

  @debug "using all factors connected to frontals and attached to separator"
  fctsyms = getFactorsAmongVariablesOnly(dfg, varlist, unused=true )
  # filter only factors connected to frontals (for upward)
  frtfcts = union(map(x->ls(dfg, x), getCliqFrontalVarIds(cliq))...)
  fctsyms = intersect(fctsyms, frtfcts)

  getCliqueData(cliq).potentials = fctsyms
  getCliqueData(cliq).partialpotential = Vector{Bool}(undef, length(fctsyms))
  fcts = map(x->getFactor(dfg, x), fctsyms)
  getCliqueData(cliq).partialpotential = map(x->isPartial(x), fcts)
  for fct in fcts
    getSolverData(fct).potentialused = true
  end

  @debug "finding all frontals for down WIP"
  ffctsyms = getCliqFactorsFromFrontals(dfg, cliq, Symbol[], inseparator=false, unused=false, solvable=solvable)
  # fnsyms = getCliqVarsWithFrontalNeighbors(csmc.dfg, csmc.cliq)
  getCliqueData(cliq).dwnPotentials = ffctsyms
  getCliqueData(cliq).dwnPartialPotential = map(x->isPartial(getFactor(dfg,x)), ffctsyms )

  nothing
end

getCliquePotentials(cliq::TreeBelief) = getCliqueData(cliq).potentials

function cliqPotentialIDs(cliq::TreeClique)
  potIDs = Symbol[]
  for idfct in getCliqueData(cliq).potentials
    push!(potIDs, idfct)
  end
  return potIDs
end

"""
    $SIGNATURES

Collect and returl all child clique separator variables.
"""
function collectSeparators(bt::AbstractBayesTree, cliq::TreeClique)::Vector{Symbol}
  allseps = Symbol[]
  for child in childCliqs(bt, cliq)#tree
      allseps = [allseps; getCliqueData(child).separatorIDs]
  end
  return allseps
end

"""
    $SIGNATURES

Return boolean matrix of factor by variable (row by column) associations within
clique, corresponds to order presented by `getCliqFactorIds` and `getCliqAllVarIds`.
"""
function getCliqAssocMat(cliq::TreeClique)
  getCliqueData(cliq).cliqAssocMat
end

"""
    $SIGNATURES

Return boolean matrix of upward message singletons (i.e. marginal priors) from
child cliques.  Variable order corresponds to `getCliqAllVarIds`.
"""
getCliqMsgMat(cliq::TreeClique) = getCliqueData(cliq).cliqMsgMat

"""
    $SIGNATURES

Return boolean matrix of factor variable associations for a clique, optionally
including (`showmsg::Bool=true`) the upward message singletons.  Variable order
corresponds to `getCliqAllVarIds`.
"""
function getCliqMat(cliq::TreeClique; showmsg::Bool=true)
  assocMat = getCliqAssocMat(cliq)
  msgMat = getCliqMsgMat(cliq)
  mat = showmsg ? [assocMat;msgMat] : assocMat
  return mat
end

"""
    $SIGNATURES

Get `cliq` separator (a.k.a. conditional) variable ids`::Symbol`.
"""
getCliqSeparatorVarIds(cliqdata::Union{BayesTreeNodeData, PackedBayesTreeNodeData})::Vector{Symbol} = cliqdata.separatorIDs
getCliqSeparatorVarIds(cliq::TreeClique)::Vector{Symbol} = getCliqSeparatorVarIds(getCliqueData(cliq))

"""
    $SIGNATURES

Get `cliq` potentials (factors) ids`::Int`.
"""
getCliqFactorIds(cliqdata::BayesTreeNodeData)::Vector{Symbol} = cliqdata.potentials
getCliqFactorIds(cliq::TreeClique)::Vector{Symbol} = getCliqFactorIds(getCliqueData(cliq))

"""
    $SIGNATURES

Get all `cliq` variable ids`::Symbol`.

Related

getCliqVarIdsAll, getCliqAllFactIds, getCliqVarsWithFrontalNeighbors
"""
function getCliqAllVarIds(cliq::TreeClique)::Vector{Symbol}
  frtl = getCliqFrontalVarIds(cliq)
  cond = getCliqSeparatorVarIds(cliq)
  union(frtl,cond)
end


"""
    $SIGNATURES

Return all variables (including frontal factor neighbors).

Dev Notes
- TODO needs to be refactored and optimized.

Related

getCliqAllVarIds
"""
function getCliqVarsWithFrontalNeighbors(fgl::G,
                                         cliq::TreeClique;
                                         solvable::Int=1  ) where {G <: AbstractDFG}
  #
  frtl = getCliqFrontalVarIds(cliq)
  cond = getCliqSeparatorVarIds(cliq)
  syms = Symbol[]
  union!(syms,Symbol.(frtl))
  union!(syms,Symbol.(cond))

  # TODO Can we trust factors are frontal connected?
  ffcs = union( map(x->ls(fgl, x, solvable=solvable), frtl)... )
  # @show ffcs = getCliqueData(cliq).potentials
  neig = union( map(x->ls(fgl, x, solvable=solvable), ffcs)... )
  union!(syms, Symbol.(neig))
  return syms
end

"""
    $SIGNATURES

Get all `cliq` variable ids`::Symbol`.

Related

getCliqAllVarIds, getCliqAllFactIds
"""
getCliqVarIdsAll(cliq::TreeClique)::Vector{Symbol} = getCliqAllVarIds(cliq::TreeClique)

"""
    $SIGNATURES

Get all `cliq` factor ids`::Symbol`.

DEPRECATED, use getCliqFactorIdsAll instead.

Related

getCliqVarIdsAll, getCliqFactors
"""
getCliqFactorIdsAll(cliqd::BayesTreeNodeData) = cliqd.potentials
getCliqFactorIdsAll(cliq::TreeClique) = getCliqFactorIdsAll(getCliqueData(cliq))
getCliqFactorIdsAll(treel::AbstractBayesTree, frtl::Symbol) = getCliqFactorIdsAll(getClique(treel, frtl))

const getCliqFactors = getCliqFactorIdsAll

"""
    $SIGNATURES

Get all `cliq` factor ids`::Symbol`.

DEPRECATED, use getCliqFactorIdsAll instead.

Related

getCliqVarIdsAll
"""
function getCliqAllFactIds(cliqd::BayesTreeNodeData)
    @warn "getCliqAllFactIds deprecated, use getCliqFactorIdsAll instead."
    return getCliqFactorIdsAll(cliqd)
end

function getCliqAllFactIds(cliq::TreeClique)
    @warn "getCliqAllFactIds deprecated, use getCliqFactorIdsAll instead."
    return getCliqFactorIdsAll(getCliqueData(cliq))
end


"""
    $SIGNATURES

Get all `cliq` variable labels as `::Symbol`.
"""
function getCliqAllVarSyms(dfg::G, cliq::TreeClique)::Vector{Symbol} where G <: AbstractDFG
  # Symbol[getSym(dfg, varid) for varid in getCliqAllVarIds(cliq)]
  @warn "deprecated getCliqAllVarSyms, use getCliqAllVarIds instead."
  getCliqAllVarIds(cliq) # not doing all frontals
end


"""
    $SIGNATURES

Return the number of factors associated with each variable in `cliq`.
"""
getCliqNumAssocFactorsPerVar(cliq::TreeClique)::Vector{Int} = sum(getCliqAssocMat(cliq), dims=1)[:]


"""
    $SIGNATURES

Get variable ids`::Int` with prior factors associated with this `cliq`.

Notes:
- does not include any singleton messages from upward or downward message passing.
"""
function getCliqVarIdsPriors(cliq::TreeClique,
                             allids::Vector{Symbol}=getCliqAllVarIds(cliq),
                             partials::Bool=true  )::Vector{Symbol}
  # get ids with prior factors associated with this cliq
  amat = getCliqAssocMat(cliq)
  prfcts = sum(amat, dims=2) .== 1

  # remove partials as per request
  !partials ? nothing : (prfcts .&= getCliqueData(cliq).partialpotential)

  # return variable ids in `mask`
  mask = sum(amat[prfcts[:],:], dims=1)[:] .> 0
  return allids[mask]
end


"""
    $SIGNATURES

Get `cliq` variable IDs with singleton factors -- i.e. both in clique priors and up messages.
"""
function getCliqVarSingletons(cliq::TreeClique,
                              allids::Vector{Symbol}=getCliqAllVarIds(cliq),
                              partials::Bool=true  )::Vector{Symbol}
  # get incoming upward messages (known singletons)
  mask = sum(getCliqMsgMat(cliq),dims=1)[:] .>= 1
  upmsgids = allids[mask]

  # get ids with prior factors associated with this cliq
  prids = getCliqVarIdsPriors(cliq, getCliqAllVarIds(cliq), partials)

  # return union of both lists
  return union(upmsgids, prids)
end


"""
    $SIGNATURES

Get each clique subgraph association matrix.
"""
function compCliqAssocMatrices!(dfg::G, bt::AbstractBayesTree, cliq::TreeClique) where G <: AbstractDFG
  frtl = getCliqFrontalVarIds(cliq)
  cond = getCliqSeparatorVarIds(cliq)
  inmsgIDs = collectSeparators(bt, cliq)
  potIDs = cliqPotentialIDs(cliq)
  # Construct associations matrix here
  # matrix has variables are columns, and messages/constraints as rows
  cols = [frtl;cond]
  getCliqueData(cliq).inmsgIDs = inmsgIDs
  getCliqueData(cliq).potIDs = potIDs
  cliqAssocMat = Array{Bool,2}(undef, length(potIDs), length(cols))
  cliqMsgMat = Array{Bool,2}(undef, length(inmsgIDs), length(cols))
  fill!(cliqAssocMat, false)
  fill!(cliqMsgMat, false)
  for j in 1:length(cols)
    for i in 1:length(inmsgIDs)
      if cols[j] == inmsgIDs[i]
        cliqMsgMat[i,j] = true
      end
    end
    for i in 1:length(potIDs)
      idfct = getCliqueData(cliq).potentials[i]
      if idfct == potIDs[i] # sanity check on clique potentials ordering
        # TODO int and symbol compare is no good
        for vertidx in getVariableOrder(DFG.getFactor(dfg, idfct)) #.fncargvID
          if vertidx == cols[j]
            cliqAssocMat[i,j] = true
          end
        end
      else
        prtslperr("compCliqAssocMatrices! -- potential ID ordering was lost")
      end
    end
  end
  getCliqueData(cliq).cliqAssocMat = cliqAssocMat
  getCliqueData(cliq).cliqMsgMat = cliqMsgMat
  nothing
end


function countSkips(bt::AbstractBayesTree)
  skps = 0
  for cliq in getCliques(bt)
    m = getCliqMat(cliq[2])
    mi = map(Int,m)
    skps += sum(map(Int,sum(mi, dims=1) .== 1))
  end
  return skps
end

function skipThroughMsgsIDs(cliq::TreeClique)
  cliqdata = getCliqueData(cliq)
  numfrtl1 = floor(Int,length(cliqdata.frontalIDs)+1)
  condAssocMat = cliqdata.cliqAssocMat[:,numfrtl1:end]
  condMsgMat = cliqdata.cliqMsgMat[:,numfrtl1:end]
  mat = [condAssocMat;condMsgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabM = sum(map(Int,condMsgMat),dims=1) .== 1
  mab = mab .& mabM
  # rang = 1:size(condMsgMat,2)
  msgidx = cliqdata.separatorIDs[vec(collect(mab))]
  return msgidx
end

function directPriorMsgIDs(cliq::TreeClique)
  frtl = getCliqueData(cliq).frontalIDs
  sepr = getCliqueData(cliq).separatorIDs
  cols = [frtl;sepr]
  mat = getCliqMat(cliq, showmsg=true)
  singr = sum(map(Int,mat),dims=2) .== 1
  rerows = collect(1:length(singr))
  b = vec(collect(singr))
  rerows2 = rerows[b]
  sumsrAc = sum(map(Int,mat[rerows2,:]),dims=1)
  sumc = sum(map(Int,mat),dims=1)
  pmSkipCols = (sumsrAc - sumc) .== 0
  return cols[vec(collect(pmSkipCols))]
end

function directFrtlMsgIDs(cliq::TreeClique)
  numfrtl = length(getCliqueData(cliq).frontalIDs)
  frntAssocMat = getCliqueData(cliq).cliqAssocMat[:,1:numfrtl]
  frtlMsgMat = getCliqueData(cliq).cliqMsgMat[:,1:numfrtl]
  mat = [frntAssocMat; frtlMsgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabM = sum(map(Int,frtlMsgMat),dims=1) .== 1
  mab = mab .& mabM
  return getCliqueData(cliq).frontalIDs[vec(collect(mab))]
end

function directAssignmentIDs(cliq::TreeClique)
  # NOTE -- old version been included in iterated variable stack
  assocMat = getCliqueData(cliq).cliqAssocMat
  msgMat = getCliqueData(cliq).cliqMsgMat
  mat = [assocMat;msgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabA = sum(map(Int,assocMat),dims=1) .== 1
  mab = mab .& mabA
  # TODO -- use proper accessor methods
  frtl = getCliqueData(cliq).frontalIDs
  sepr = getCliqueData(cliq).separatorIDs
  cols = [frtl;sepr]
  return cols[vec(collect(mab))]
  # also calculate how which are conditionals
end

function mcmcIterationIDs(cliq::TreeClique)
  @debug "mcmcIterationIDs\n" cliq.id getCliqFrontalVarIds(cliq) getCliqSeparatorVarIds(cliq)
  mat = getCliqMat(cliq)
  @debug "getCliqMat" mat
  # assocMat = getCliqueData(cliq).cliqAssocMat
  # msgMat = getCliqueData(cliq).cliqMsgMat
  # mat = [assocMat;msgMat];
  sum(sum(map(Int,mat),dims=1)) == 0 ? error("mcmcIterationIDs -- unaccounted variables") : nothing
  mab = 1 .< sum(map(Int,mat),dims=1)
  cols = getCliqAllVarIds(cliq)

  # must also include "direct variables" connected through projection only
  directvars = directAssignmentIDs(cliq)
  usset = union(directvars, cols[vec(collect(mab))])
  # NOTE -- fix direct vs itervar issue, DirectVarIDs against Iters should also Iter
  # NOTE -- using direct then mcmcIter ordering to prioritize non-msg vars first
  return setdiff(usset, getCliqueData(cliq).directPriorMsgIDs)
end

function getCliqMatVarIdx(cliq::TreeClique, varid::Symbol, allids=getCliqAllVarIds(cliq) )
  len = length(allids)
  [1:len;][allids .== varid][1]
end

"""
    $SIGNATURES

Determine and return order list of variable ids required for minibatch Gibbs iteration inside `cliq`.

Notes
- Singleton factors (priors and up messages) back of the list
- least number of associated factor variables earlier in list
- Same as getCliqVarSolveOrderUp
"""
function mcmcIterationIdsOrdered(cliq::TreeClique)
  # get unordered iter list
  alliter = mcmcIterationIDs(cliq)

  # get all singletons
  allsings = getCliqVarSingletons(cliq)
  singletonvars = intersect(alliter, allsings)

  # get all non-singleton iters
  nonsinglvars = setdiff(alliter, singletonvars)

  # sort nonsingletons ascending number of factors
  mat = getCliqMat(cliq)
  lenfcts = sum(mat, dims=1)
  nonslen = zeros(length(nonsinglvars))
  for i in 1:length(nonsinglvars)
    varid = nonsinglvars[i]
    varidx = getCliqMatVarIdx(cliq, varid)
    nonslen[i] = lenfcts[varidx]
  end
  p = sortperm(nonslen)
  ascnons = nonsinglvars[p]

  # sort singleton vars ascending number of factors
  singslen = zeros(length(singletonvars))
  for i in 1:length(singletonvars)
    varid = singletonvars[i]
    varidx = getCliqMatVarIdx(cliq, varid)
    singslen[i] = lenfcts[varidx]
  end
  p = sortperm(singslen)
  ascsing = singletonvars[p]

  return [ascnons; ascsing]
end

"""
    $SIGNATURES

Determine and return order list of variable ids required for minibatch Gibbs iteration inside `cliq`.

Notes
- Singleton factors (priors and up messages) back of the list
- least number of associated factor variables earlier in list
- Same as mcmcIterationIdsOrdered
"""
function getCliqVarSolveOrderUp(cliq::TreeClique)
  return mcmcIterationIdsOrdered(cliq)
end

"""
    $(SIGNATURES)

Prepare the variable IDs for nested clique Gibbs mini-batch calculations, by assembing these clique data fields:
- `directPriorMsgIDs`
- `directvarIDs`
- `itervarIDs`
- `msgskipIDs`
- `directFrtlMsgIDs`

"""
function setCliqMCIDs!(cliq::TreeClique)
  getCliqueData(cliq).directPriorMsgIDs = directPriorMsgIDs(cliq)

  # NOTE -- directvarIDs are combined into itervarIDs
  getCliqueData(cliq).directvarIDs = directAssignmentIDs(cliq)
  # TODO find itervarIDs that have upward child singleton messages and update them last in iter list
  getCliqueData(cliq).itervarIDs = mcmcIterationIdsOrdered(cliq)

  getCliqueData(cliq).msgskipIDs = skipThroughMsgsIDs(cliq)
  getCliqueData(cliq).directFrtlMsgIDs = directFrtlMsgIDs(cliq)

  # TODO add initialization sequence var id list too

  nothing
end


# post order tree traversal and build potential functions
function buildCliquePotentials(dfg::G, bt::AbstractBayesTree, cliq::TreeClique; solvable::Int=1) where G <: AbstractDFG
    for child in childCliqs(bt, cliq)#tree
        buildCliquePotentials(dfg, bt, child)
    end
    @debug "Get potentials $(getLabel(cliq))"
    setCliqPotentials!(dfg, bt, cliq, solvable=solvable)

    compCliqAssocMatrices!(dfg, bt, cliq)
    setCliqMCIDs!(cliq)

    nothing
end

"""
    $(SIGNATURES)

Return a vector of child cliques to `cliq`.
"""
function childCliqs(treel::MetaBayesTree, cliq::TreeClique)
  cliqKey = treel.bt[:cliqId][cliq.id]
  childcliqs = TreeClique[]
  for cIdx in MetaGraphs.outneighbors(treel.bt, cliqKey)
    push!(childcliqs, get_prop(treel.bt, cIdx, :clique))
  end
  return childcliqs
end

"""
    $(SIGNATURES)

Return a vector of child cliques to `cliq`.
"""
getChildren(treel::AbstractBayesTree, frtsym::Symbol) = childCliqs(treel, frtsym)
getChildren(treel::AbstractBayesTree, cliq::TreeClique) = childCliqs(treel, cliq)

"""
    $SIGNATURES
Get edges to children cliques
"""
function getEdgesChildren(tree::MetaBayesTree, cliqkey::Int)
  [MetaGraphs.Edge(cliqkey, chkey) for chkey in MetaGraphs.outneighbors(tree.bt, cliqkey)]
end

getEdgesChildren(tree::MetaBayesTree, cliq::TreeClique) = getEdgesChildren(tree, tree.bt[:cliqId][cliq.id])

"""
    $SIGNATURES
Get edges to parent clique
"""
function getEdgesParent(tree::MetaBayesTree, cliqkey::Int)
  [MetaGraphs.Edge(pkey, cliqkey) for pkey in MetaGraphs.inneighbors(tree.bt, cliqkey)]
end

getEdgesParent(tree::MetaBayesTree, cliq::TreeClique) = getEdgesParent(tree, tree.bt[:cliqId][cliq.id])

"""
    $SIGNATURES

Return a vector of all siblings to a clique, which defaults to not `inclusive` the calling `cliq`.
"""
function getCliqSiblings(treel::AbstractBayesTree, cliq::TreeClique, inclusive::Bool=false)::Vector{TreeClique}
  prnt = getParent(treel, cliq)
  if length(prnt) > 0
    allch = getChildren(treel, prnt[1])
  end
  if inclusive
    return allch
  end
  sibs = TreeClique[]
  for ch in allch
    if ch.id != cliq.id
      push!(sibs, ch)
    end
  end
  return sibs
end

"""
    $(SIGNATURES)

Return `cliq`'s parent clique.
"""
function parentCliq(treel::MetaBayesTree, cliq::TreeClique)
  cliqKey = treel.bt[:cliqId][cliq.id]
  parentcliqs = TreeClique[]
  for pIdx in  MetaGraphs.inneighbors(treel.bt, cliqKey)
    push!(parentcliqs, get_prop(treel.bt, pIdx, :clique))
  end
  return parentcliqs
end

"""
    $(SIGNATURES)

Return number of cliques in a tree.
"""
getNumCliqs(tree::MetaBayesTree) = MetaGraphs.nv(tree.bt)


"""
    $(SIGNATURES)

Return `cliq`'s parent clique.
"""
getParent(treel::AbstractBayesTree, afrontal::Union{Symbol, TreeClique}) = parentCliq(treel, afrontal)

"""
    $SIGNATURES

Return one symbol (a frontal variable) from each clique in the `::BayesTree`.

Notes
- Frontal variables only occur once in a clique per tree, therefore is a unique identifier.

Related:

whichCliq, printCliqHistorySummary
"""
function getTreeAllFrontalSyms(::AbstractDFG, tree::AbstractBayesTree)
  cliqs = getCliques(tree)
  syms = Vector{Symbol}(undef, length(cliqs))
  for (id,cliq) in cliqs
    syms[id] = getCliqFrontalVarIds(cliq)[1]
  end
  return syms
end


# import DistributedFactorGraphs: getVariableOrder
#
# getVariableOrder(treel::AbstractBayesTree)::Vector{Symbol} = treel.variableOrder

"""
    $SIGNATURES

Return the variable elimination order stored in a tree object.
"""
getEliminationOrder(treel::AbstractBayesTree) = treel.variableOrder


"""
    $SIGNATURES

EXPERIMENTAL, Save a Bayes (Junction) tree object to file.

Notes
- Converts and saves to JLD2 format a set of `PackedBayesTreeNodeData` objects.
- IIF issue #481

Related

IIF.loadTree, DFG.saveDFG, DFG.loadDFG, JLD2.@save, JLD2.@load
"""
function saveTree(treel::AbstractBayesTree,
                  filepath=joinpath("/tmp","caesar","savetree.jld2") )
  #
  savetree = deepcopy(treel)
  for i in 1:length(getCliques(savetree))
    if  getCliqueData(savetree, i) isa BayesTreeNodeData
      setCliqueData!(getClique(savetree, i), convert(PackedBayesTreeNodeData, getCliqueData(savetree, i)))
    end
  end

  JLD2.@save filepath savetree
  return filepath
end

function saveTree(treeArr::Vector{T},
                  filepath=joinpath("/tmp","caesar","savetrees.jld2") ) where T <: AbstractBayesTree
  #
  savetree = deepcopy(treeArr)
  for savtre in savetree, i in 1:length(getCliques(savtre))
    if getCliqueData(savtre, i) isa BayesTreeNodeData
      setCliqueData!(getClique(savtre,i), convert(PackedBayesTreeNodeData, getCliqueData(savtre, i)))
    end
  end

  JLD2.@save filepath savetree
  return filepath
end

"""
    $SIGNATURES

EXPERIMENTAL, Save a Bayes (Junction) tree object to file.

Notes
- Converts and saves to JLD2 format a set of `PackedBayesTreeNodeData` objects.
- IIF issue #481

Related

IIF.saveTree, DFG.saveDFG, DFG.loadDFG, JLD2.@save, JLD2.@load
"""
function loadTree(filepath=joinpath("/tmp","caesar","savetree.jld2"))
  data = @load filepath savetree

  # convert back to a type that which could not be serialized by JLD2
  if savetree isa Vector
      for savtre in savetree, i in 1:length(getCliques(savtre))
        if getCliqueData(savtre, i) isa PackedBayesTreeNodeData
          setCliqueData!(getClique(savtre, i), convert(BayesTreeNodeData, getCliqueData(savtre, i)))
        end
      end
  else
    for i in 1:length(getCliques(savetree))
      if getCliqueData(savetree, i) isa PackedBayesTreeNodeData
        setCliqueData!(getClique(savetree, i), convert(BayesTreeNodeData, getCliqueData(savetree, i)))
      end
    end
  end

  # return loaded and converted tree
  return savetree
end

"""
    $SIGNATURES

Return Tuple of number cliques (Marginalized, Reused).
"""
function calcCliquesRecycled(tree::AbstractBayesTree)
  numMarg = 0
  numReused = 0
  numBoth = 0

  for (key, cliq) in tree.cliques
    numReused += getCliqueData(cliq).isCliqReused ? 1 : 0
    numMarg += getCliqueData(cliq).allmarginalized ? 1 : 0
    numBoth += getCliqueData(cliq).allmarginalized && getCliqueData(cliq).isCliqReused ? 1 : 0
  end

  return length(tree.cliques), numMarg, numReused, numBoth
end


## Tree Reuse

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
    clq = TreeClique(-1)
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

