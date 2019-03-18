
"""
$(TYPEDEF)
"""
mutable struct BayesTreeNodeData
  frontalIDs::Vector{Int}
  conditIDs::Vector{Int}
  inmsgIDs::Vector{Int}
  potIDs::Vector{Int} # this is likely redundant TODO -- remove
  potentials::Vector{Int}
  cliqAssocMat::Array{Bool,2}
  cliqMsgMat::Array{Bool,2}
  directvarIDs::Vector{Int}
  directFrtlMsgIDs::Vector{Int}
  msgskipIDs::Vector{Int}
  itervarIDs::Vector{Int}
  directPriorMsgIDs::Vector{Int}
  debug
  debugDwn
  upMsg::Dict{Symbol, BallTreeDensity}
  dwnMsg::Dict{Symbol, BallTreeDensity}
  BayesTreeNodeData() = new()
  BayesTreeNodeData(x...) = new(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16])
end

# TODO -- this should be a constructor
function emptyBTNodeData()
  BayesTreeNodeData(Int[],Int[],Int[],
                    Int[],Int[],Array{Bool}(undef, 0,0),
                    Array{Bool}(undef, 0,0),Int[],Int[],
                    Int[],Int[],Int[],
                    nothing, nothing,
                    Dict{Symbol, BallTreeDensity}(:null => AMP.manikde!(zeros(1,1), [1.0;], (:Euclid,))),
                    Dict{Symbol, BallTreeDensity}(:null => AMP.manikde!(zeros(1,1), [1.0;], (:Euclid,))) )
end

# BayesTree declarations
"""
$(TYPEDEF)
"""
mutable struct BayesTree
  bt
  btid::Int
  cliques::Dict{Int,Graphs.ExVertex}
  frontals::Dict{String,Int}
end

function emptyBayesTree()
    bt =   BayesTree(Graphs.inclist(Graphs.ExVertex,is_directed=true),
                     0,
                     Dict{Int,Graphs.ExVertex}(),
                     #[],
                     Dict{AbstractString, Int}())
    return bt
end

# create a new clique
function addClique!(bt::BayesTree, fg::FactorGraph, varID::Int, condIDs::Array{Int}=Int[])
  bt.btid += 1
  clq = Graphs.add_vertex!(bt.bt, ExVertex(bt.btid,string("Clique",bt.btid)))
  bt.cliques[bt.btid] = clq

  clq.attributes["label"] = ""
  # Specific data container
  clq.attributes["data"] = emptyBTNodeData()

  appendClique!(bt, bt.btid, fg, varID, condIDs)
  return clq
end

# generate the label for particular clique -- graphviz drawing
function makeCliqueLabel(fgl::FactorGraph, bt::BayesTree, clqID::Int)
  clq = bt.cliques[clqID]
  flbl = ""
  clbl = ""
  for fr in clq.attributes["data"].frontalIDs
    flbl = string(flbl,localapi.getvertex(fgl,fr).attributes["label"], ",") #fgl.v[fr].
  end
  for cond in clq.attributes["data"].conditIDs
    clbl = string(clbl, localapi.getvertex(fgl,cond).attributes["label"], ",") # fgl.v[cond].
  end
  clq.attributes["label"] = string(flbl, ": ", clbl)
end

# add a conditional ID to clique
function appendConditional(bt::BayesTree, clqID::Int, condIDs::Array{Int,1})
  clq = bt.cliques[clqID]
  clq.attributes["data"].conditIDs = union(clq.attributes["data"].conditIDs, condIDs)
end

# Add a new frontal variable to clique
function appendClique!(bt::BayesTree, clqID::Int, fg::FactorGraph, varID::Int, condIDs::Array{Int,1}=Int[])
  clq = bt.cliques[clqID]
  var = localapi.getvertex(fg, varID) # fg.v[varID]
  # add frontal variable
  push!(clq.attributes["data"].frontalIDs,varID)
  # total dictionary of frontals for easy access
  bt.frontals[var.attributes["label"]] = clqID#bt.btid

  appendConditional(bt, clqID, condIDs)
  makeCliqueLabel(fg, bt, clqID)
  nothing
end


# instantiate a new child clique in the tree
function newChildClique!(bt::BayesTree, fg::FactorGraph, CpID::Int, varID::Int, Sepj::Array{Int,1})
  chclq = addClique!(bt, fg, varID, Sepj)
  parent = bt.cliques[CpID]
  # Staying with Graphs.jl for tree in first stage
  edge = Graphs.make_edge(bt.bt, parent, chclq)
  Graphs.add_edge!(bt.bt, edge)

  return chclq
end

# post order tree traversal and build potential functions
function findCliqueFromFrontal(bt::BayesTree, frtlID::Int)
  for cliqPair in bt.cliques
    id = cliqPair[1]
    cliq = cliqPair[2]
    for frtl in cliq.attributes["data"].frontalIDs
      if frtl == frtlID
        return cliq
      end
    end
  end
  error("Clique with desired frontal ID not found")
end


# eliminate a variable for new
function newPotential(tree::BayesTree, fg::FactorGraph, var::Int, prevVar::Int, p::Array{Int,1})
    firvert = localapi.getvertex(fg,var)
    if (length(getData(firvert).separator) == 0)
      if (length(tree.cliques) == 0)
        addClique!(tree, fg, var)
      else
        appendClique!(tree, 1, fg, var) # add to root
      end
    else
      Sj = getData(firvert).separator
      # find parent clique Cp that containts the first eliminated variable of Sj as frontal
      firstelim = 99999999999
      for s in Sj
        temp = something(findfirst(isequal(s), p), 0) # findfirst(p, s)
        if (temp < firstelim)
          firstelim = temp
        end
      end
      felbl = localapi.getvertex(fg, p[firstelim]).attributes["label"]
      CpID = tree.frontals[felbl]
      # look to add this conditional to the tree
      unFC = union(tree.cliques[CpID].attributes["data"].frontalIDs, tree.cliques[CpID].attributes["data"].conditIDs)
      if (sort(unFC) == sort(Sj))
        appendClique!(tree, CpID, fg, var)
      else
        newChildClique!(tree, fg, CpID, var, Sj)
      end
    end
end

# build the whole tree in batch format
function buildTree!(tree::BayesTree, fg::FactorGraph, p::Array{Int,1})
  rp = reverse(p,dims=1) # flipdim(p, 1)
  prevVar = 0
  for var in rp
    newPotential(tree, fg, var, prevVar, p)
    prevVar = var
  end
end

function showTree(;filepath::String="/tmp/bt.pdf",
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

function drawTree(treel::BayesTree;
                  show::Bool=false,                  # must remain false for stability and automated use in solver
                  filepath::String="/tmp/bt.pdf",
                  viewerapp::String="evince"  )
  #
  fext = split(filepath, '.')[end]
  fpwoext = split(filepath, '.')[end-1]
  fid = IOStream("")
  try
    fid = open("$(fpwoext).dot","w+")
    write(fid,to_dot(treel.bt))
    close(fid)
    run(`dot $(fpwoext).dot -T$(fext) -o $(filepath)`)
  catch ex
    @warn ex
    @show stacktrace()
  finally
    close(fid)
  end

  show ? showTree(viewerapp=viewerapp, filepath=filepath) : nothing
end



## Find batch belief propagation solution
function prepBatchTree!(fg::FactorGraph;
                        ordering::Symbol=:qr,
                        drawpdf::Bool=false,
                        show::Bool=false,
                        filepath::String="/tmp/bt.pdf",
                        viewerapp::String="evince"  )
  #
  p = IncrementalInference.getEliminationOrder(fg, ordering=ordering)
  println()
  fge = deepcopy(fg)
  println("Building Bayes net...")
  buildBayesNet!(fge, p)

  tree = emptyBayesTree()
  buildTree!(tree, fge, p)

  println("Bayes Net")
  # sleep(0.1)
  #fid = open("bn.dot","w+")
  #write(fid,to_dot(fge.bn))
  #close(fid)

  # Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
  println("Bayes Tree")
  if drawpdf
    drawTree(tree, show=show, filepath=filepath, viewerapp=viewerapp)
  end

  # GraphViz.Graph(to_dot(tree.bt))
  #Michael reference 3sig -- x2l1x4x3    x1|x2

  println("Find potential functions for each clique")
  cliq = tree.cliques[1] # start at the root
  buildCliquePotentials(fg, tree, cliq); # fg does not have the marginals as fge does

  # now update all factor graph vertices used for this tree
  for v in vertices(fg.g)
    dlapi.updatevertex!(fg, v)
  end

  return tree
end

function resetData!(vdata::VariableNodeData)
  vdata.eliminated = false
  vdata.BayesNetOutVertIDs = Int[]
  vdata.BayesNetVertID = 0
  vdata.separator = Int[]
  nothing
end

function resetData!(vdata::FunctionNodeData)
  vdata.eliminated = false
  vdata.potentialused = false
  nothing
end

function resetFactorGraphNewTree!(fg::FactorGraph)
  for v in vertices(fg.g)
    resetData!(getData(v))
    localapi.updatevertex!(fg, v)
  end

  nothing
end

"""
    $(SIGNATURES)

Build a completely new Bayes (Junction) tree, after first wiping clean all temporary state in fg from a possibly pre-existing tree.
"""
function wipeBuildNewTree!(fg::FactorGraph;
                           ordering=:qr,
                           drawpdf=false,
                           show::Bool=false,
                           filepath::String="/tmp/bt.pdf",
                           viewerapp::String="evince"  )
  #
  resetFactorGraphNewTree!(fg);
  return prepBatchTree!(fg, ordering=ordering, drawpdf=drawpdf, show=show, filepath=filepath, viewerapp=viewerapp);
end

"""
    $(SIGNATURES)

Return the Graphs.ExVertex node object that represents a clique in the Bayes (Junction) tree, as defined by one of the frontal variables `frt`.
"""
function whichCliq(bt::BayesTree, frt::T) where {T <: AbstractString}
    bt.cliques[bt.frontals[frt]]
end
whichCliq(bt::BayesTree, frt::Symbol) = whichCliq(bt, string(frt))

"""
    $SIGNATURES

Return the Graphs.ExVertex node object that represents a clique in the Bayes (Junction) tree, as defined by one of the frontal variables `frt`.
"""
getCliq(bt::BayesTree, frt::Symbol) = whichCliq(bt, string(frt))



"""
    $(SIGNATURES)

Set the upward passing message for Bayes (Junction) tree clique `cliql`.
"""
function setUpMsg!(cliql::ExVertex, msgs::Dict{Symbol, BallTreeDensity})
  getData(cliql).upMsg = msgs
end

"""
    $(SIGNATURES)

Set the downward passing message for Bayes (Junction) tree clique `cliql`.
"""
function setDwnMsg!(cliql::ExVertex, msgs::Dict{Symbol, BallTreeDensity})
  getData(cliql).dwnMsg = msgs
end

"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
function upMsg(cliq::Graphs.ExVertex)
  getData(cliq).upMsg
end
function upMsg(btl::BayesTree, sym::Symbol)
  upMsg(whichCliq(btl, sym))
end

"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
getUpMsgs(btl::BayesTree, sym::Symbol) = upMsg(btl, sym)
getUpMsgs(cliql::Graphs.ExVertex) = upMsg(cliql)



"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
function dwnMsg(cliq::Graphs.ExVertex)
  getData(cliq).dwnMsg
end
function dwnMsg(btl::BayesTree, sym::Symbol)
  upMsg(whichCliq(btl, sym))
end

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
getDwnMsgs(btl::BayesTree, sym::Symbol) = dwnMsg(btl, sym)
getDwnMsgs(cliql::Graphs.ExVertex) = dwnMsg(cliql)


function appendUseFcts!(usefcts, lblid::Int, fct::Graphs.ExVertex, fid::Int)
  for tp in usefcts
    if tp == fct.index
      return nothing
    end
  end
  tpl = fct.index
  push!(usefcts, tpl )
  nothing
end



function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
    frtl = cliq.attributes["data"].frontalIDs
    cond = cliq.attributes["data"].conditIDs
    allids = [frtl;cond]
    alldimIDs = Int[]
    for fid in frtl
      alldimIDs = [alldimIDs; getData(localapi.getvertex(fg,fid)).dimIDs]
    end
    for cid in cond
      alldimIDs = [alldimIDs; getData(localapi.getvertex(fg,cid)).dimIDs]
    end

    for fid in frtl
        usefcts = []
        for fct in localapi.outneighbors(fg, localapi.getvertex(fg,fid))
            if getData(fct).potentialused!=true
                loutn = localapi.outneighbors(fg, fct)
                if length(loutn)==1
                    appendUseFcts!(usefcts, fg.IDs[Symbol(loutn[1].label)], fct, fid)
                    # TODO -- make update vertex call
                    fct.attributes["data"].potentialused = true
                    localapi.updatevertex!(fg, fct)
                end
                for sepSearch in loutn
                    sslbl = Symbol(sepSearch.label)
                    if (fg.IDs[sslbl] == fid)
                        continue # skip the fid itself
                    end
                    sea = findmin(abs.(allids .- fg.IDs[sslbl]))
                    if sea[1]==0.0
                        appendUseFcts!(usefcts, fg.IDs[sslbl], fct, fid)
                        # usefcts = [usefcts;(fg.IDs[sslbl], fct, fid)]
                        fct.attributes["data"].potentialused = true #fct.attributes["potentialused"] = true
                        localapi.updatevertex!(fg, fct)
                    end
                end
            end
        end
        cliq.attributes["data"].potentials=union(cliq.attributes["data"].potentials,usefcts)
    end
    return nothing
end

function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, chkcliq::Int)
    getCliquePotentials!(fg, bt.cliques[chkcliq])
end

function cliqPotentialIDs(cliq::Graphs.ExVertex)
  potIDs = Int[]
  for idfct in cliq.attributes["data"].potentials
    push!(potIDs,idfct)
  end
  return potIDs
end

function collectSeparators(bt::BayesTree, cliq::Graphs.ExVertex)
  allseps = Int[]
  for child in out_neighbors(cliq, bt.bt)#tree
      allseps = [allseps; child.attributes["data"].conditIDs]
  end
  return allseps
end

function compCliqAssocMatrices!(fgl::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
  frtl = cliq.attributes["data"].frontalIDs
  cond = cliq.attributes["data"].conditIDs
  inmsgIDs = collectSeparators(bt, cliq)
  potIDs = cliqPotentialIDs(cliq)
  # Construct associations matrix here
  # matrix has variables are columns, and messages/constraints as rows
  cols = [frtl;cond]
  cliq.attributes["data"].inmsgIDs = inmsgIDs
  cliq.attributes["data"].potIDs = potIDs
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
      idfct = cliq.attributes["data"].potentials[i]
      if idfct == potIDs[i] # sanity check on clique potentials ordering
        for vertidx in getData(getVert(fgl, idfct, api=localapi)).fncargvID
        # for vertidx in getData(getVertNode(fgl, idfct)).fncargvID
          if vertidx == cols[j]
            cliqAssocMat[i,j] = true
          end
        end
      else
        prtslperr("compCliqAssocMatrices! -- potential ID ordering was lost")
      end
    end
  end
  cliq.attributes["data"].cliqAssocMat = cliqAssocMat
  cliq.attributes["data"].cliqMsgMat = cliqMsgMat
  nothing
end

function getCliqAssocMat(cliq::Graphs.ExVertex)
  cliq.attributes["data"].cliqAssocMat
end
function getCliqMsgMat(cliq::Graphs.ExVertex)
  cliq.attributes["data"].cliqMsgMat
end
function getCliqMat(cliq::Graphs.ExVertex; showmsg=true)
  assocMat = getCliqAssocMat(cliq)
  msgMat = getCliqMsgMat(cliq)
  mat = showmsg ? [assocMat;msgMat] : assocMat
  return mat
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
  msgidx = cliqdata.conditIDs[vec(collect(mab))]
  return msgidx
end

function directPriorMsgIDs(cliq::Graphs.ExVertex)
  frtl = getData(cliq).frontalIDs
  cond = getData(cliq).conditIDs
  cols = [frtl;cond]
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
  frtl = getData(cliq).frontalIDs
  cond = getData(cliq).conditIDs
  cols = [frtl;cond]
  return cols[vec(collect(mab))]
  # also calculate how which are conditionals
end

function mcmcIterationIDs(cliq::Graphs.ExVertex)
  assocMat = getData(cliq).cliqAssocMat
  msgMat = getData(cliq).cliqMsgMat
  mat = [assocMat;msgMat];
  sum(sum(map(Int,mat),dims=1)) == 0 ? error("mcmcIterationIDs -- unaccounted variables") : nothing
  mab = 1 .< sum(map(Int,mat),dims=1)
  frtl = getData(cliq).frontalIDs
  cond = getData(cliq).conditIDs
  cols = [frtl;cond]

  return cols[vec(collect(mab))]
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

  # NOTE -- combined larger iter group
  getData(cliq).directvarIDs = directAssignmentIDs(cliq)
  # getData(cliq).itervarIDs = mcmcIterationIDs(cliq)
  # NOTE -- fix direct vs itervar issue, DirectVarIDs against Iters should also Iter
  # NOTE -- using direct then mcmcIter ordering to prioritize non-msg vars first
  usset = union(directAssignmentIDs(cliq), mcmcIterationIDs(cliq))
  getData(cliq).itervarIDs = setdiff(usset, getData(cliq).directPriorMsgIDs)
  getData(cliq).msgskipIDs = skipThroughMsgsIDs(cliq)
  getData(cliq).directFrtlMsgIDs = directFrtlMsgIDs(cliq)

  # TODO find itervarIDs that have upward child singleton messages and update them last in iter list


  nothing
end


# post order tree traversal and build potential functions
function buildCliquePotentials(fg::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
    for child in out_neighbors(cliq, bt.bt)#tree
        buildCliquePotentials(fg, bt, child)
    end
    @info "Get potentials $(cliq.attributes["label"])"
    getCliquePotentials!(fg, bt, cliq);

    compCliqAssocMatrices!(fg, bt, cliq);
    setCliqMCIDs!(cliq);

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
getParent(treel::BayesTree, afrontal::Symbol) = parentCliq(treel, afrontal)
