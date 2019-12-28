export getVariableOrder

"""
    $SIGNATURES

Get the frontal variable IDs `::Int` for a given clique in a Bayes (Junction) tree.
"""
getCliqFrontalVarIds(cliqdata::Union{BayesTreeNodeData, PackedBayesTreeNodeData})::Vector{Symbol} = cliqdata.frontalIDs
getCliqFrontalVarIds(cliq::Graphs.ExVertex)::Vector{Symbol} = getCliqFrontalVarIds(getData(cliq))

"""
    $SIGNATURES

Get the frontal variable IDs `::Int` for a given clique in a Bayes (Junction) tree.
"""
getFrontals(cliqd::Union{Graphs.ExVertex,BayesTreeNodeData})::Vector{Symbol} = getCliqFrontalVarIds(cliqd)


"""
    $SIGNATURES

Create a new clique.
"""
function addClique!(bt::BayesTree, dfg::G, varID::Symbol, condIDs::Array{Symbol}=Symbol[])::ExVertex where G <: AbstractDFG
  bt.btid += 1
  clq = Graphs.add_vertex!(bt.bt, ExVertex(bt.btid,string("Clique",bt.btid)))
  bt.cliques[bt.btid] = clq

  clq.attributes["label"] = ""

  # Specific data container
  setData!(clq, emptyBTNodeData())
  # clq.attributes["data"] = emptyBTNodeData()

  appendClique!(bt, bt.btid, dfg, varID, condIDs)
  return clq
end

"""
    $SIGNATURES

Generate the label for particular clique (used by graphviz for visualization).
"""
function makeCliqueLabel(dfg::G, bt::BayesTree, clqID::Int)::String where G <: AbstractDFG
  clq = bt.cliques[clqID]
  flbl = ""
  clbl = ""
  for fr in getData(clq).frontalIDs
    flbl = string(flbl, DFG.getVariable(dfg,fr).label, ",")
  end
  for sepr in getData(clq).separatorIDs
    clbl = string(clbl, DFG.getVariable(dfg,sepr).label, ",")
  end
  clq.attributes["label"] = string(flbl, ": ", clbl)
end

"""
    $SIGNATURES

Add the separator for the newly created clique.
"""
function appendSeparatorToClique(bt::BayesTree, clqID::Int, seprIDs::Array{Symbol,1})
  union!(getData(bt.cliques[clqID]).separatorIDs, seprIDs)
  nothing
end

"""
    $SIGNATURES

Add a new frontal variable to clique.

DevNotes
- TODO, define what "conditionals" are CLEARLY!!
"""
function appendClique!(bt::BayesTree,
                       clqID::Int,
                       dfg::AbstractDFG,
                       varID::Symbol,
                       seprIDs::Array{Symbol,1}=Symbol[] )::Nothing
  #
  clq = bt.cliques[clqID]
  var = DFG.getVariable(dfg, varID)

  # add frontal variable
  push!(getData(clq).frontalIDs, varID)

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
function newChildClique!(bt::BayesTree, dfg::AbstractDFG, CpID::Int, varID::Symbol, Sepj::Array{Symbol,1})
  # physically create the new clique
  chclq = addClique!(bt, dfg, varID, Sepj)
  parent = bt.cliques[CpID]
  # Staying with Graphs.jl for tree in first stage
  edge = Graphs.make_edge(bt.bt, parent, chclq)
  Graphs.add_edge!(bt.bt, edge)

  return chclq
end

"""
    $SIGNATURES

Post order tree traversal and build potential functions
"""
function findCliqueFromFrontal(bt::BayesTree, frtlID::Int)
  for cliqPair in bt.cliques
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
                                          Sj=solverData(firvert).separator)::DFGVariable
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
function newPotential(tree::BayesTree, dfg::G, var::Symbol, elimorder::Array{Symbol,1}) where G <: AbstractDFG
    firvert = DFG.getVariable(dfg,var)
    # no parent
    if (length(solverData(firvert).separator) == 0)
      # if (length(tree.cliques) == 0)
        # create new root
        addClique!(tree, dfg, var)
      # else
      #   # add to root
      #   @warn "root append clique is happening"
      #   appendClique!(tree, 1, dfg, var)
      # end
    else
      # find parent clique Cp that containts the first eliminated variable of Sj as frontal
      Sj = solverData(firvert).separator
      felbl = identifyFirstEliminatedSeparator(dfg, elimorder, firvert, Sj).label
      # get clique id of first eliminated frontal
      CpID = tree.frontals[felbl]
      # look to add this conditional to the tree
      cliq = tree.cliques[CpID]
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
function buildTree!(tree::BayesTree, dfg::AbstractDFG, elimorder::Array{Symbol,1})
  revorder = reverse(elimorder,dims=1) # flipdim(p, 1), fixing #499
  # prevVar = 0
  for var in revorder
    @info "Adding $var to tree..."
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
function showTree(;filepath::String="/tmp/caesar/bt.pdf",
                   viewerapp::String="evince"  )
  #
  try
    @async run(`$(viewerapp) $(filepath)`)
  catch ex
    @warn "not able to show via $(viewerapp) $(filepath)"
    @show ex
    @show stacktrace()
  end
end

"""
    $SIGNATURES

Draw the Bayes (Junction) tree by means of `.dot` files and `pdf` reader.

Notes
- Uses system install of graphviz.org.
- Can also use Linux tool `xdot`.
"""
function drawTree(treel::BayesTree;
                  show::Bool=false,                  # must remain false for stability and automated use in solver
                  filepath::String="/tmp/caesar/bt.pdf",
                  viewerapp::String="evince",
                  imgs::Bool=false )
  #
  fext = filepath[(end-2):end] #split(filepath, '.')[end]
  fpwoext = filepath[1:(end-4)]# split(filepath, '.')[end-1]
  mkpath(dirname(fpwoext))

  # modify a deepcopy
  btc = deepcopy(treel)
  for (cid, cliq) in btc.cliques
    if imgs
      firstlabel = split(cliq.attributes["label"],',')[1]
      spyCliqMat(cliq, suppressprint=true) |> exportimg("/tmp/$firstlabel.png")
      cliq.attributes["image"] = "/tmp/$firstlabel.png"
      cliq.attributes["label"] = ""
    end
    delete!(cliq.attributes, "data")
  end

  fid = IOStream("")
  try
    fid = open("$(fpwoext).dot","w+")
    write(fid,to_dot(btc.bt))
    close(fid)
    run(`dot $(fpwoext).dot -T $(fext) -o $(filepath)`)
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

Draw the Bayes (junction) tree with LaTeX labels by means of `.dot` and `.tex`
files.

Notes
- Uses system install of graphviz.org.
- Uses external python `dot2tex` tool (`pip install dot2tex`).

Related:

drawTree
"""
function generateTexTree(treel::BayesTree;
                         filepath::String="/tmp/caesar/bayes/bt")
    btc = deepcopy(treel)
    for (cid, cliq) in btc.cliques
        label = cliq.attributes["label"]

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
        cliq.attributes["label"] = newlabel
    end

    # Use new labels to produce `.dot` and `.tex` files.
    fid = IOStream("")
    try
        mkpath(joinpath((split(filepath,'/')[1:(end-1)])...) )
        fid = open("$(filepath).dot","w+")
        write(fid, to_dot(btc.bt))
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
function buildTreeFromOrdering!(dfg::G,
                                p::Vector{Symbol};
                                drawbayesnet::Bool=false,
                                maxparallel::Int=50,
                                solvable::Int=1 ) where G <: InMemoryDFGTypes
  #
  println()
  fge = deepcopy(dfg)
  @info "Building Bayes net..."
  buildBayesNet!(fge, p, maxparallel=maxparallel, solvable=solvable)

  @info "Staring the Bayes tree construction from Bayes net"
  tree = emptyBayesTree()
  tree.variableOrder = p
  buildTree!(tree, fge, p)

  if drawbayesnet
    println("Bayes Net")
    sleep(0.1)
    fid = open("bn.dot","w+")
    write(fid,to_dot(fge.bn))
    close(fid)
  end

  println("Find potential functions for each clique")
  cliq = tree.cliques[1] # start at the root
  buildCliquePotentials(dfg, tree, cliq, solvable=solvable); # fg does not have the marginals as fge does

  return tree
end


"""
    $SIGNATURES

Build Bayes/Junction/Elimination tree from a given variable ordering.
"""
function buildTreeFromOrdering!(dfg::DFG.AbstractDFG,
                                p::Vector{Symbol};
                                drawbayesnet::Bool=false,
                                maxparallel::Int=50  )
  #
  println()

  @info "Copying to a local DFG"
  fge = InMemDFGType(params=getSolverParams(dfg))#GraphsDFG{SolverParams}(params=SolverParams())
  DistributedFactorGraphs._copyIntoGraph!(dfg, fge, union(getVariableIds(dfg), getFactorIds(dfg)), true)

  println("Building Bayes net from cloud...")
  buildBayesNet!(fge, p, maxparallel=maxparallel)

  tree = emptyBayesTree()
  tree.variableOrder = p
  buildTree!(tree, fge, p)

  if drawbayesnet
    println("Bayes Net")
    sleep(0.1)
    fid = open("bn.dot","w+")
    write(fid,to_dot(fge.bn))
    close(fid)
  end

  println("Find potential functions for each clique")
  cliq = tree.cliques[1] # start at the root
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
                        ordering::Symbol=:qr,
                        variableOrder::Union{Nothing, Vector{Symbol}}=nothing,
                        drawpdf::Bool=false,
                        show::Bool=false,
                        filepath::String="/tmp/caesar/bt.pdf",
                        viewerapp::String="evince",
                        imgs::Bool=false,
                        drawbayesnet::Bool=false,
                        maxparallel::Int=50 )
  #
  p = variableOrder != nothing ? variableOrder : getEliminationOrder(dfg, ordering=ordering)

  # for debuggin , its useful to have the variable ordering
  if drawpdf
    ispath(getLogPath(dfg)) ? nothing : Base.mkpath(getLogPath(dfg))
    open(joinLogPath(dfg,"variableOrder.txt"), "a") do io
      writedlm(io, string.(reshape(p,1,:)), ',')
    end
  end

  tree = buildTreeFromOrdering!(dfg, p, drawbayesnet=drawbayesnet, maxparallel=maxparallel)

  # GraphViz.Graph(to_dot(tree.bt))
  # Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
  #Michael reference 3sig -- x2l1x4x3    x1|x2
  println("Bayes Tree")
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
function resetFactorGraphNewTree!(dfg::AbstractDFG)::Nothing
  for v in DFG.getVariables(dfg)
    resetData!(solverData(v))
  end
  for f in DFG.getFactors(dfg)
    resetData!(solverData(f))
  end
  nothing
end

"""
    $SIGNATURES

Reset factor graph and build a new tree from the provided variable ordering `p`.
"""
function resetBuildTreeFromOrder!(fgl::AbstractDFG, p::Vector{Symbol})::BayesTree
  resetFactorGraphNewTree!(fgl)
  return buildTreeFromOrdering!(fgl, p)
end

"""
    $(SIGNATURES)

Build a completely new Bayes (Junction) tree, after first wiping clean all
temporary state in fg from a possibly pre-existing tree.

Related:

buildTreeFromOrdering!
"""
function wipeBuildNewTree!(dfg::G;
                           ordering::Symbol=:qr,
                           drawpdf::Bool=false,
                           show::Bool=false,
                           filepath::String="/tmp/caesar/bt.pdf",
                           viewerapp::String="evince",
                           imgs::Bool=false,
                           maxparallel::Int=50,
                           variableOrder::Union{Nothing, Vector{Symbol}}=nothing  )::BayesTree where G <: AbstractDFG
  #
  resetFactorGraphNewTree!(dfg);
  return prepBatchTree!(dfg, variableOrder=variableOrder, ordering=ordering, drawpdf=drawpdf, show=show, filepath=filepath, viewerapp=viewerapp, imgs=imgs, maxparallel=maxparallel);
end

"""
    $(SIGNATURES)

Return the Graphs.ExVertex node object that represents a clique in the Bayes
(Junction) tree, as defined by one of the frontal variables `frt<:AbstractString`.

Notes
- Frontal variables only occur once in a clique per tree, therefore is a unique identifier.

Related:

getCliq, getTreeAllFrontalSyms
"""
getCliq(bt::BayesTree, frt::Symbol) = bt.cliques[bt.frontals[frt]]

"""
    $(SIGNATURES)

Return the Graphs.ExVertex node object that represents a clique in the Bayes
(Junction) tree, as defined by one of the frontal variables `frt<:AbstractString`.

Notes
- Frontal variables only occur once in a clique per tree, therefore is a unique identifier.

Related:

getCliq, getTreeAllFrontalSyms
"""
whichCliq(bt::BayesTree, frt::Symbol) = getCliq(bt, frt)
whichCliq(bt::BayesTree, frt::T) where {T <: AbstractString} = whichCliq(bt, Symbol(string(frt)))


"""
    $SIGNATURES

Return boolean on whether the frontal variable `frt::Symbol` exists somewhere in the `::BayesTree`.
"""
hasCliq(bt::BayesTree, frt::Symbol)::Bool = haskey(bt.frontals, frt)

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
getCliqDepth(tree::BayesTree, sym::Symbol)::Int = getCliqDepth(tree, getCliq(tree, sym))



"""
    $(SIGNATURES)

Set the upward passing message for Bayes (Junction) tree clique `cliql`.

Dev Notes
- TODO setUpMsg! should also set inferred dimension
"""
function setUpMsg!(cliql::ExVertex, msgs::Dict{Symbol, BallTreeDensity})
  @error "setUpMsg!, use inferred dimension version instead"
  getData(cliql).upMsg = msgs
end

function setUpMsg!(cliql::ExVertex, msgs::TempBeliefMsg) #Dict{Symbol, Tuple{BallTreeDensity, Float64}}
  # ms = Dict{Symbol, BallTreeDensity}()
  # for (id, val) in msgs
  #   ms[id] = val[1]
  # end
  getData(cliql).upMsg = msgs #ms
  nothing
end

"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
getUpMsgs(cliql::Graphs.ExVertex) = getData(cliql).upMsg
getUpMsgs(btl::BayesTree, sym::Symbol) = getUpMsgs(whichCliq(btl, sym))

"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
getCliqMsgsUp(cliql::Graphs.ExVertex) = upMsg(cliql)
getCliqMsgsUp(treel::BayesTree, frt::Symbol) = getCliqMsgsUp(getCliq(treel, frt))

"""
    $(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""
function setDwnMsg!(cliql::ExVertex, msgs::TempBeliefMsg) #Dict{Symbol, BallTreeDensity}
  getData(cliql).dwnMsg = msgs
end

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
getDwnMsgs(cliql::Graphs.ExVertex) = getData(cliql).dwnMsg
getDwnMsgs(btl::BayesTree, sym::Symbol) = getDwnMsgs(whichCliq(btl, sym))

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
getCliqMsgsDown(cliql::Graphs.ExVertex) = getDwnMsgs(cliql)


function appendUseFcts!(usefcts,
                        lblid::Symbol,
                        fct::DFGFactor )
                        # fid::Symbol )
  #
  union!(usefcts, Symbol(fct.label))
  # for tp in usefcts
  #   if tp == fct.index
  #     return nothing
  #   end
  # end
  # tpl = fct.label
  # push!(usefcts, tpl )
  nothing
end

"""
    $SIGNATURES

Return list of factors which depend only on variables in variable list in factor
graph -- i.e. among variables.

Notes
-----
* `unused::Bool=true` will disregard factors already used -- i.e. disregard where `potentialused=true`
"""
function getFactorsAmongVariablesOnly(dfg::G,
                                      varlist::Vector{Symbol};
                                      unused::Bool=true  ) where G <: AbstractDFG
  # collect all factors attached to variables
  prefcts = Symbol[]
  for var in varlist
    union!(prefcts, DFG.ls(dfg, var))
  end

  almostfcts = Symbol[]
  if unused
    # now check if those factors have already been added
    for fct in prefcts
      vert = DFG.getFactor(dfg, fct)
      if !solverData(vert).potentialused
        push!(almostfcts, fct)
      end
    end
  else
    almostfcts = prefcts
  end

  # Select factors that have all variables in this clique var list
  usefcts = Symbol[]
  for fct in almostfcts
    if length(setdiff(DFG.getNeighbors(dfg, fct), varlist)) == 0
      push!(usefcts, fct)
    end
  end

  return usefcts
end


"""
    $SIGNATURES

Get all factors connected to frontal variables

Dev Notes
- Why not just do this, `ffcs = union( map(x->ls(fgl, x), frtl) )`?
"""
function getCliqFactorsFromFrontals(fgl::G,
                                    cliq::Graphs.ExVertex,
                                    varlist::Vector{Symbol};
                                    inseparator::Bool=true,
                                    unused::Bool=true,
                                    solvable::Int=1  ) where {G <: AbstractDFG}
  #
  frtls = getData(cliq).frontalIDs
  seprs = getData(cliq).separatorIDs
  allids = [frtls;seprs]

  usefcts = Symbol[]
  for frsym in frtls
    # usefcts = Int[]
    for fctid in ls(fgl, frsym)
        fct = getFactor(fgl, fctid)
        if !unused || !solverData(fct).potentialused
            loutn = ls(fgl, fctid, solvable=solvable)
            # deal with unary factors
            if length(loutn)==1
                union!(usefcts, Symbol[Symbol(fct.label);])
                # appendUseFcts!(usefcts, loutn, fct) # , frsym)
                solverData(fct).potentialused = true
            end
            # deal with n-ary factors
            for sep in loutn
                if sep == frsym
                    continue # skip the frsym itself
                end
                insep = sep in allids
                if !inseparator || insep
                    union!(usefcts, Symbol[Symbol(fct.label);])
                    solverData(fct).potentialused = true
                    if !insep
                        @info "cliq=$(cliq.index) adding factor that is not in separator, $sep"
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
function isPartial(fct::DFGFactor)  #fct::Graphs.ExVertex
  fcf = solverData(fct).fnc.usrfnc!
  isPartial(fcf)
end

"""
    $SIGNATURES

Determine and set the potentials for a particular `cliq` in the Bayes (Junction) tree.
"""
function setCliqPotentials!(dfg::G,
                            bt::BayesTree,
                            cliq::Graphs.ExVertex;
                            solvable::Int=1  ) where G <: AbstractDFG
  #
  varlist = getCliqVarIdsAll(cliq)

  @info "using all factors connected to frontals and attached to separator"
  fctsyms = getFactorsAmongVariablesOnly(dfg, varlist, unused=true )
  # filter only factors connected to frontals (for upward)
  frtfcts = union(map(x->ls(dfg, x), getCliqFrontalVarIds(cliq))...)
  fctsyms = intersect(fctsyms, frtfcts)

  getData(cliq).potentials = fctsyms
  getData(cliq).partialpotential = Vector{Bool}(undef, length(fctsyms))
  fcts = map(x->getFactor(dfg, x), fctsyms)
  getData(cliq).partialpotential = map(x->isPartial(x), fcts)
  for fct in fcts
    solverData(fct).potentialused = true
  end

  @info "finding all frontals for down WIP"
  ffctsyms = getCliqFactorsFromFrontals(dfg, cliq, Symbol[], inseparator=false, unused=false, solvable=solvable)
  # fnsyms = getCliqVarsWithFrontalNeighbors(csmc.dfg, csmc.cliq)
  getData(cliq).dwnPotentials = ffctsyms
  getData(cliq).dwnPartialPotential = map(x->isPartial(getFactor(dfg,x)), ffctsyms )

  nothing
end

function getCliquePotentials!(dfg::G,
                            bt::BayesTree,
                            cliq::Graphs.ExVertex  ) where G <: AbstractDFG
  #
  @warn "getCliquePotentials! deprecated, use getCliqPotentials instead."
  getCliqPotentials(dfg, bt, cliq)
end

function cliqPotentialIDs(cliq::Graphs.ExVertex)
  potIDs = Symbol[]
  for idfct in getData(cliq).potentials
    push!(potIDs, idfct)
  end
  return potIDs
end

"""
    $SIGNATURES

Collect and returl all child clique separator variables.
"""
function collectSeparators(bt::BayesTree, cliq::Graphs.ExVertex)::Vector{Symbol}
  allseps = Symbol[]
  for child in out_neighbors(cliq, bt.bt)#tree
      allseps = [allseps; getData(child).separatorIDs]
  end
  return allseps
end

"""
    $SIGNATURES

Return boolean matrix of factor by variable (row by column) associations within
clique, corresponds to order presented by `getCliqFactorIds` and `getCliqAllVarIds`.
"""
function getCliqAssocMat(cliq::Graphs.ExVertex)
  getData(cliq).cliqAssocMat
end

"""
    $SIGNATURES

Return boolean matrix of upward message singletons (i.e. marginal priors) from
child cliques.  Variable order corresponds to `getCliqAllVarIds`.
"""
getCliqMsgMat(cliq::Graphs.ExVertex) = getData(cliq).cliqMsgMat

"""
    $SIGNATURES

Return boolean matrix of factor variable associations for a clique, optionally
including (`showmsg::Bool=true`) the upward message singletons.  Variable order
corresponds to `getCliqAllVarIds`.
"""
function getCliqMat(cliq::Graphs.ExVertex; showmsg::Bool=true)
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
getCliqSeparatorVarIds(cliq::Graphs.ExVertex)::Vector{Symbol} = getCliqSeparatorVarIds(getData(cliq))

"""
    $SIGNATURES

Get `cliq` potentials (factors) ids`::Int`.
"""
getCliqFactorIds(cliqdata::BayesTreeNodeData)::Vector{Symbol} = cliqdata.potentials
getCliqFactorIds(cliq::Graphs.ExVertex)::Vector{Symbol} = getCliqFactorIds(getData(cliq))

"""
    $SIGNATURES

Get all `cliq` variable ids`::Symbol`.

Related

getCliqVarIdsAll, getCliqAllFactIds, getCliqVarsWithFrontalNeighbors
"""
function getCliqAllVarIds(cliq::Graphs.ExVertex)::Vector{Symbol}
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
                                         cliq::Graphs.ExVertex;
                                         solvable::Int=1  ) where {G <: AbstractDFG}
  #
  frtl = getCliqFrontalVarIds(cliq)
  cond = getCliqSeparatorVarIds(cliq)
  syms = Symbol[]
  union!(syms,Symbol.(frtl))
  union!(syms,Symbol.(cond))

  # TODO Can we trust factors are frontal connected?
  ffcs = union( map(x->ls(fgl, x, solvable=solvable), frtl)... )
  # @show ffcs = getData(cliq).potentials
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
getCliqVarIdsAll(cliq::Graphs.ExVertex)::Vector{Symbol} = getCliqAllVarIds(cliq::Graphs.ExVertex)

"""
    $SIGNATURES

Get all `cliq` factor ids`::Symbol`.

DEPRECATED, use getCliqFactorIdsAll instead.

Related

getCliqVarIdsAll, getCliqFactors
"""
getCliqFactorIdsAll(cliqd::BayesTreeNodeData) = cliqd.potentials
getCliqFactorIdsAll(cliq::Graphs.ExVertex) = getCliqFactorIdsAll(getData(cliq))
getCliqFactorIdsAll(treel::BayesTree, frtl::Symbol) = getCliqFactorIdsAll(getCliq(treel, frtl))

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

function getCliqAllFactIds(cliq::Graphs.ExVertex)
    @warn "getCliqAllFactIds deprecated, use getCliqFactorIdsAll instead."
    return getCliqFactorIdsAll(getData(cliq))
end


"""
    $SIGNATURES

Get all `cliq` variable labels as `::Symbol`.
"""
function getCliqAllVarSyms(dfg::G, cliq::Graphs.ExVertex)::Vector{Symbol} where G <: AbstractDFG
  # Symbol[getSym(dfg, varid) for varid in getCliqAllVarIds(cliq)]
  @warn "deprecated getCliqAllVarSyms, use getCliqAllVarIds instead."
  getCliqAllVarIds(cliq) # not doing all frontals
end

"""
    $SIGNATURES

Return dictionary of down messages consisting of all frontal and separator beliefs of this clique.

Notes:
- Fetches numerical results from `subdfg` as dictated in `cliq`.
"""
function getCliqDownMsgsAfterDownSolve(subdfg::G, cliq::Graphs.ExVertex)::TempBeliefMsg where G <: AbstractDFG
  # Dict{Symbol, BallTreeDensity}
  # where the return msgs are contained
  container = TempBeliefMsg() # Dict{Symbol,BallTreeDensity}()

  # go through all msgs one by one
  for sym in getCliqAllVarIds(cliq)
    container[sym] = (getKDE(subdfg, sym), getVariableInferredDim(subdfg, sym))
  end

  # return the result
  return container
end



"""
    $SIGNATURES

Get variable ids`::Int` with prior factors associated with this `cliq`.

Notes:
- does not include any singleton messages from upward or downward message passing.
"""
function getCliqVarIdsPriors(cliq::Graphs.ExVertex,
                             allids::Vector{Symbol}=getCliqAllVarIds(cliq),
                             partials::Bool=true  )::Vector{Symbol}
  # get ids with prior factors associated with this cliq
  amat = getCliqAssocMat(cliq)
  prfcts = sum(amat, dims=2) .== 1

  # remove partials as per request
  !partials ? nothing : (prfcts .&= getData(cliq).partialpotential)

  # return variable ids in `mask`
  mask = sum(amat[prfcts[:],:], dims=1)[:] .> 0
  return allids[mask]
end

"""
    $SIGNATURES

Return the number of factors associated with each variable in `cliq`.
"""
getCliqNumAssocFactorsPerVar(cliq::Graphs.ExVertex)::Vector{Int} = sum(getCliqAssocMat(cliq), dims=1)[:]



"""
    $SIGNATURES

Get `cliq` variable IDs with singleton factors -- i.e. both in clique priors and up messages.
"""
function getCliqVarSingletons(cliq::Graphs.ExVertex,
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
function compCliqAssocMatrices!(dfg::G, bt::BayesTree, cliq::Graphs.ExVertex) where G <: AbstractDFG
  frtl = getCliqFrontalVarIds(cliq)
  cond = getCliqSeparatorVarIds(cliq)
  inmsgIDs = collectSeparators(bt, cliq)
  potIDs = cliqPotentialIDs(cliq)
  # Construct associations matrix here
  # matrix has variables are columns, and messages/constraints as rows
  cols = [frtl;cond]
  getData(cliq).inmsgIDs = inmsgIDs
  getData(cliq).potIDs = potIDs
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
      idfct = getData(cliq).potentials[i]
      if idfct == potIDs[i] # sanity check on clique potentials ordering
        # TODO int and symbol compare is no good
        for vertidx in solverData(DFG.getFactor(dfg, idfct)).fncargvID
          if vertidx == cols[j]
            cliqAssocMat[i,j] = true
          end
        end
      else
        prtslperr("compCliqAssocMatrices! -- potential ID ordering was lost")
      end
    end
  end
  getData(cliq).cliqAssocMat = cliqAssocMat
  getData(cliq).cliqMsgMat = cliqMsgMat
  nothing
end


function countSkips(bt::BayesTree)
  skps = 0
  for cliq in bt.cliques
    m = getCliqMat(cliq[2])
    mi = map(Int,m)
    skps += sum(map(Int,sum(mi, dims=1) .== 1))
  end
  return skps
end

function skipThroughMsgsIDs(cliq::Graphs.ExVertex)
  cliqdata = getData(cliq)
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

function directPriorMsgIDs(cliq::Graphs.ExVertex)
  frtl = getData(cliq).frontalIDs
  sepr = getData(cliq).separatorIDs
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

function directFrtlMsgIDs(cliq::Graphs.ExVertex)
  numfrtl = length(getData(cliq).frontalIDs)
  frntAssocMat = getData(cliq).cliqAssocMat[:,1:numfrtl]
  frtlMsgMat = getData(cliq).cliqMsgMat[:,1:numfrtl]
  mat = [frntAssocMat; frtlMsgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabM = sum(map(Int,frtlMsgMat),dims=1) .== 1
  mab = mab .& mabM
  return getData(cliq).frontalIDs[vec(collect(mab))]
end

function directAssignmentIDs(cliq::Graphs.ExVertex)
  # NOTE -- old version been included in iterated variable stack
  assocMat = getData(cliq).cliqAssocMat
  msgMat = getData(cliq).cliqMsgMat
  mat = [assocMat;msgMat];
  mab = sum(map(Int,mat),dims=1) .== 1
  mabA = sum(map(Int,assocMat),dims=1) .== 1
  mab = mab .& mabA
  # TODO -- use proper accessor methods
  frtl = getData(cliq).frontalIDs
  sepr = getData(cliq).separatorIDs
  cols = [frtl;sepr]
  return cols[vec(collect(mab))]
  # also calculate how which are conditionals
end

function mcmcIterationIDs(cliq::Graphs.ExVertex)
  mat = getCliqMat(cliq)
  # assocMat = getData(cliq).cliqAssocMat
  # msgMat = getData(cliq).cliqMsgMat
  # mat = [assocMat;msgMat];

  sum(sum(map(Int,mat),dims=1)) == 0 ? error("mcmcIterationIDs -- unaccounted variables") : nothing
  mab = 1 .< sum(map(Int,mat),dims=1)
  cols = getCliqAllVarIds(cliq)

  # must also include "direct variables" connected through projection only
  directvars = directAssignmentIDs(cliq)
  usset = union(directvars, cols[vec(collect(mab))])
  # NOTE -- fix direct vs itervar issue, DirectVarIDs against Iters should also Iter
  # NOTE -- using direct then mcmcIter ordering to prioritize non-msg vars first
  return setdiff(usset, getData(cliq).directPriorMsgIDs)
end

function getCliqMatVarIdx(cliq::Graphs.ExVertex, varid::Symbol, allids=getCliqAllVarIds(cliq) )
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
function mcmcIterationIdsOrdered(cliq::Graphs.ExVertex)
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
function getCliqVarSolveOrderUp(cliq::Graphs.ExVertex)
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
function setCliqMCIDs!(cliq::Graphs.ExVertex)
  getData(cliq).directPriorMsgIDs = directPriorMsgIDs(cliq)

  # NOTE -- directvarIDs are combined into itervarIDs
  getData(cliq).directvarIDs = directAssignmentIDs(cliq)
  # TODO find itervarIDs that have upward child singleton messages and update them last in iter list
  getData(cliq).itervarIDs = mcmcIterationIdsOrdered(cliq)  #mcmcIterationIDs(cliq)

  getData(cliq).msgskipIDs = skipThroughMsgsIDs(cliq)
  getData(cliq).directFrtlMsgIDs = directFrtlMsgIDs(cliq)

  # TODO add initialization sequence var id list too

  nothing
end


# post order tree traversal and build potential functions
function buildCliquePotentials(dfg::G, bt::BayesTree, cliq::Graphs.ExVertex; solvable::Int=1) where G <: AbstractDFG
    for child in out_neighbors(cliq, bt.bt)#tree
        buildCliquePotentials(dfg, bt, child)
    end
    @info "Get potentials $(cliq.attributes["label"])"
    setCliqPotentials!(dfg, bt, cliq, solvable=solvable)

    compCliqAssocMatrices!(dfg, bt, cliq)
    setCliqMCIDs!(cliq)

    nothing
end

"""
    $(SIGNATURES)

Return a vector of child cliques to `cliq`.
"""
function childCliqs(treel::BayesTree, cliq::Graphs.ExVertex)
    childcliqs = Vector{Graphs.ExVertex}(undef, 0)
    for cl in Graphs.out_neighbors(cliq, treel.bt)
        push!(childcliqs, cl)
    end
    return childcliqs
end
function childCliqs(treel::BayesTree, frtsym::Symbol)
  childCliqs(treel,  whichCliq(treel, frtsym))
end

"""
    $(SIGNATURES)

Return a vector of child cliques to `cliq`.
"""
getChildren(treel::BayesTree, frtsym::Symbol) = childCliqs(treel, frtsym)
getChildren(treel::BayesTree, cliq::Graphs.ExVertex) = childCliqs(treel, cliq)

"""
    $SIGNATURES

Return a vector of all siblings to a clique, which defaults to not `inclusive` the calling `cliq`.
"""
function getCliqSiblings(treel::BayesTree, cliq::Graphs.ExVertex, inclusive::Bool=false)::Vector{Graphs.ExVertex}
  prnt = getParent(treel, cliq)
  if length(prnt) > 0
    allch = getChildren(treel, prnt[1])
  end
  if inclusive
    return allch
  end
  sibs = Graphs.ExVertex[]
  for ch in allch
    if ch.index != cliq.index
      push!(sibs, ch)
    end
  end
  return sibs
end

"""
    $(SIGNATURES)

Return `cliq`'s parent clique.
"""
function parentCliq(treel::BayesTree, cliq::Graphs.ExVertex)
    Graphs.in_neighbors(cliq, treel.bt)
end
function parentCliq(treel::BayesTree, frtsym::Symbol)
  parentCliq(treel,  whichCliq(treel, frtsym))
end

"""
    $(SIGNATURES)

Return `cliq`'s parent clique.
"""
getParent(treel::BayesTree, afrontal::Union{Symbol, Graphs.ExVertex}) = parentCliq(treel, afrontal)

"""
    $SIGNATURES

Return one symbol (a frontal variable) from each clique in the `::BayesTree`.

Notes
- Frontal variables only occur once in a clique per tree, therefore is a unique identifier.

Related:

whichCliq, printCliqHistorySummary
"""
function getTreeAllFrontalSyms(fgl::G, tree::BayesTree) where G <: AbstractDFG
  cliqs = tree.cliques
  syms = Vector{Symbol}(undef, length(cliqs))
  for (id,cliq) in cliqs
    syms[id] = getCliqFrontalVarIds(cliq)[1]
  end
  return syms
end

"""
    $SIGNATURES

Get the `::Condition` variable for a clique, likely used for delaying state transitions in
state machine solver.
"""
getSolveCondition(cliq::Graphs.ExVertex) = getData(cliq).solveCondition


"""
    $SIGNATURES

Return dictionary of all up belief messages currently in a Bayes `tree`.

Related

getUpMsgs
"""
function getTreeCliqUpMsgsAll(tree::BayesTree)::Dict{Int,TempBeliefMsg}
  allUpMsgs = Dict{Int,TempBeliefMsg}()
  for (idx,cliq) in tree.cliques
    msgs = getUpMsgs(cliq)
    allUpMsgs[cliq.index] = TempBeliefMsg()
    for (lbl,msg) in msgs
      # TODO capture the inferred dimension as part of the upward propagation
      allUpMsgs[cliq.index][lbl] = msg
    end
  end
  return allUpMsgs
end

"""
    $SIGNATURES

Convert tree up messages dictionary to a new dictionary relative to variables specific messages and their depth in the tree

Notes
- Return data in `TempUpMsgPlotting` format:
    Dict{Symbol,   -- is for variable label
     Vector{       -- multiple msgs for the same variable
      Symbol,      -- Clique index
      Int,         -- Depth in tree
      BTD          -- Belief estimate
      inferredDim  -- Information count
     }
"""
function stackCliqUpMsgsByVariable(tree::BayesTree,
                                   tmpmsgs::Dict{Int, TempBeliefMsg}  )::TempUpMsgPlotting
  #
  # start of the return data structure
  stack = TempUpMsgPlotting()

  # look at all the clique level data
  for (cidx,tmpmsg) in tmpmsgs
    # look at all variables up msg from each clique
    for (sym,msgdim) in tmpmsg
      # create a new object for a particular variable if hasnt been seen before
      if !haskey(stack,sym)
        stack[sym] = Vector{Tuple{Symbol, Int, BallTreeDensity, Float64}}()
      end
      # assemble metadata
      cliq = tree.cliques[cidx]
      frt = getCliqFrontalVarIds(cliq)[1]
      # add this belief msg and meta data to vector of variable entry
      push!(stack[sym], (frt, getCliqDepth(tree, cliq),msgdim[1], msgdim[2]))
    end
  end

  return stack
end

"""
    $SIGNATURES

Return the variable order stored in a tree object.
"""
getVariableOrder(treel::BayesTree)::Vector{Symbol} = treel.variableOrder

getEliminationOrder(treel::BayesTree) = treel.variableOrder


"""
    $SIGNATURES

EXPERIMENTAL, Save a Bayes (Junction) tree object to file.

Notes
- Converts and saves to JLD2 format a set of `PackedBayesTreeNodeData` objects.
- IIF issue #481

Related

IIF.loadTree, DFG.saveDFG, DFG.loadDFG, JLD2.@save, JLD2.@load
"""
function saveTree(treel::BayesTree,
                  filepath=joinpath("/tmp","caesar","savetree.jld2") )
  #
  savetree = deepcopy(treel)
  for i in 1:length(savetree.cliques)
    if  savetree.cliques[i].attributes["data"] isa BayesTreeNodeData
      savetree.cliques[i].attributes["data"] = convert(PackedBayesTreeNodeData, savetree.cliques[i].attributes["data"])
    end
  end

  JLD2.@save filepath savetree
  return filepath
end

function saveTree(treeArr::Vector{BayesTree},
                  filepath=joinpath("/tmp","caesar","savetrees.jld2") )
  #
  savetree = deepcopy(treeArr)
  for savtre in savetree, i in 1:length(savtre.cliques)
    if savtre.cliques[i].attributes["data"] isa BayesTreeNodeData
      savtre.cliques[i].attributes["data"] = convert(PackedBayesTreeNodeData, savtre.cliques[i].attributes["data"])
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
      for savtre in savetree, i in 1:length(savtre.cliques)
        if savtre.cliques[i].attributes["data"] isa PackedBayesTreeNodeData
          savtre.cliques[i].attributes["data"] = convert(BayesTreeNodeData, savtre.cliques[i].attributes["data"])
        end
      end
  else
    for i in 1:length(savetree.cliques)
      if savetree.cliques[i].attributes["data"] isa PackedBayesTreeNodeData
        savetree.cliques[i].attributes["data"] = convert(BayesTreeNodeData, savetree.cliques[i].attributes["data"])
      end
    end
  end

  # return loaded and converted tree
  return savetree
end
