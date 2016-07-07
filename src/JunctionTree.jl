


# BayesTree declarations
type BayesTree
  bt
  btid::Int
  cliques::Dict{Int,Graphs.ExVertex}
  #edges
  frontals::Dict{ASCIIString,Int64}
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
  clq.attributes["frontalIDs"] = Int[]
  clq.attributes["conditIDs"] = Int[]

  clq.attributes["label"] = ""
  clq.attributes["potentials"] = []

  appendClique!(bt, bt.btid, fg, varID, condIDs)
  return clq
end

# generate the label for particular clique -- graphviz drawing
function makeCliqueLabel(fgl::FactorGraph, bt::BayesTree, clqID::Int)
  clq = bt.cliques[clqID]
  flbl = ""
  clbl = ""
  for fr in clq.attributes["frontalIDs"]
    flbl = string(flbl,dlapi.getvertex(fgl,fr).attributes["label"], ",") #fgl.v[fr].
  end
  for cond in clq.attributes["conditIDs"]
    clbl = string(clbl, dlapi.getvertex(fgl,cond).attributes["label"], ",") # fgl.v[cond].
  end
  clq.attributes["label"] = string(flbl, ": ", clbl)
end

# add a conditional ID to clique
function appendConditional(bt::BayesTree, clqID::Int, condIDs::Array{Int,1})
  #@show "adding conditionals", condIDs
  clq = bt.cliques[clqID]
  clq.attributes["conditIDs"] = union(clq.attributes["conditIDs"], condIDs)
end

# Add a new frontal variable to clique
function appendClique!(bt::BayesTree, clqID::Int, fg::FactorGraph, varID::Int, condIDs::Array{Int,1}=Int[])
  clq = bt.cliques[clqID]
  var = dlapi.getvertex(fg, varID) # fg.v[varID]
  # add frontal variable
  push!(clq.attributes["frontalIDs"],varID)
  # total dictionary of frontals for easy access
  bt.frontals[var.attributes["label"]] = clqID#bt.btid

  appendConditional(bt, clqID, condIDs)
  makeCliqueLabel(fg, bt, clqID)
  Union{}
end


# instantiate a new child clique in the tree
function newChildClique!(bt::BayesTree, fg::FactorGraph, CpID::Int, varID::Int, Sepj::Array{Int,1})
  chclq = addClique!(bt, fg, varID, Sepj)
  parent = bt.cliques[CpID]
  # Staying with Graphs.jl for tree in first stage
  edge = Graphs.make_edge(bt.bt, parent, chclq)
  Graphs.add_edge!(bt.bt, edge)
  # addEdge!(bt.bt, parent, chclq)

  return chclq
end

# post order tree traversal and build potential functions
function findCliqueFromFrontal(bt::BayesTree, frtlID::Int)
  for cliqPair in bt.cliques
    id = cliqPair[1]
    cliq = cliqPair[2]
    for frtl in cliq.attributes["frontalIDs"]
      if frtl == frtlID
        return cliq
      end
    end
  end
  error("Clique with desired frontal ID not found")
end


# eliminate a variable for new
function newPotential(tree::BayesTree, fg::FactorGraph, var::Int, prevVar::Int, p::Array{Int,1})
    #@show fg.v[var].attributes["label"]
    if (length(dlapi.getvertex(fg,var).attributes["data"].separator) == 0) #fg.v[var].
      if (length(tree.cliques) == 0)
        addClique!(tree, fg, var)
      else
        appendClique!(tree, 1, fg, var) # add to root
      end
    else
      Sj = dlapi.getvertex(fg,var).attributes["data"].separator  # fg.v[var].
      # Sj = fg.v[var].attributes["separator"]
      # find parent clique Cp that containts the first eliminated variable of Sj as frontal
      firstelim = 99999999999
      for s in Sj
        temp = findfirst(p, s)
        if (temp < firstelim)
          firstelim = temp
        end
      end
      felbl = dlapi.getvertex(fg,p[firstelim]).attributes["label"] # fg.v[p[firstelim]].
      CpID = tree.frontals[felbl]
      # look to add this conditional to the tree
      unFC = union(tree.cliques[CpID].attributes["frontalIDs"], tree.cliques[CpID].attributes["conditIDs"])
      if (sort(unFC) == sort(Sj))
        appendClique!(tree, CpID, fg, var)
      else
        newChildClique!(tree, fg, CpID, var, Sj)
      end
    end
end

# build the whole tree in batch format
function buildTree!(tree::BayesTree, fg::FactorGraph, p::Array{Int,1})
  rp = flipdim(p,1)
  prevVar = 0
  for var in rp
    newPotential(tree, fg, var, prevVar, p)
    prevVar = var
  end
end


marg(x...) = -1, [0.0] # TODO -- this should be removed


## Find batch belief propagation solution
function prepBatchTree!(fg::FactorGraph; ordering::Symbol=:qr,drawpdf::Bool=false)
  p = getEliminationOrder(fg, ordering=ordering)
  #@show p=[1, 3, 8, 10, 5]
  # for v in p
  #     print("$(fg.v[v].label),")
  # end
  println()

  fge = deepcopy(fg)
  buildBayesNet!(fge, p)

  tree = emptyBayesTree()
  buildTree!(tree, fge, p)

  println("Bayes Net")
  sleep(0.1)
  #fid = open("bn.dot","w+")
  #write(fid,to_dot(fge.bn))
  #close(fid)

  #GraphViz.Graph(to_dot(fge.bn))
  # Michael reference -- x2->x1, x2->x3, x2->x4, x2->l1, x4->x3, l1->x3, l1->x4
  println("Bayes Tree")
  if drawpdf
    fid = open("bt.dot","w+")
    write(fid,to_dot(tree.bt))
    close(fid)
    run(`dot bt.dot -Tpdf -o bt.pdf`)
  end

  # GraphViz.Graph(to_dot(tree.bt))
  #Michael reference 3sig -- x2l1x4x3    x1|x2

  println("Find potential functions for each clique")
  cliq = tree.cliques[1] # start at the root
  buildCliquePotentials(fg, tree, cliq);

  return tree
end

function resetFactorGraphNewTree!(fg::FactorGraph)
  for v in fg.v
    v[2].attributes["data"].eliminated = false
    # v[2].attributes["eliminated"] = false
    v[2].attributes["data"].BayesNetOutVertIDs = Int64[]
    v[2].attributes["data"].BayesNetVertID = 0
    v[2].attributes["data"].separator = Int64[]
    # v[2].attributes["BayesNetVert"] = Union{}
    # v[2].attributes["separator"] = Int64[]
    dlapi.updatevertex!(fg, v)
  end
  for f in fg.f
    f[2].attributes["data"].eliminated = false #f[2].attributes["eliminated"] = false
    f[2].attributes["data"].potentialused = false #f[2].attributes["potentialused"] = false
    dlapi.updatevertex!(fg, v)
  end

  nothing
end

function wipeBuildNewTree!(fg::FactorGraph; ordering=:qr)
  resetFactorGraphNewTree!(fg);
  return prepBatchTree!(fg, ordering=ordering);
end

function whichCliq(bt::BayesTree, frt::ASCIIString)
    bt.cliques[bt.frontals[frt]]
end


function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
    frtl = cliq.attributes["frontalIDs"]
    cond = cliq.attributes["conditIDs"]
    allids = [frtl;cond]
    alldimIDs = Int[]
    for fid in frtl
      alldimIDs = [alldimIDs;dlapi.getvertex(fg,fid).attributes["data"].dimIDs] # ;fg.v[fid].
      # alldimIDs = [alldimIDs;fg.v[fid].attributes["dimIDs"]]
    end
    for cid in cond
      alldimIDs = [alldimIDs;dlapi.getvertex(fg,cid).attributes["data"].dimIDs] # ;fg.v[cid].
      # alldimIDs = [alldimIDs;fg.v[cid].attributes["dimIDs"]]
    end

    for fid in frtl
        usefcts = []
        for fct in dlapi.outneighbors(fg, dlapi.getvertex(fg,fid)) #out_neighbors(dlapi.getvertex(fg,fid),fg.g) # (fg.v[fid],
            #println("Checking fct=$(fct.label)")
            if fct.attributes["data"].potentialused!=true #fct.attributes["potentialused"]!=true ## USED TO HAVE == AND continue end here
                loutn = dlapi.outneighbors(fg, fct) #out_neighbors(fct, fg.g)
                if length(loutn)==1
                    appendUseFcts!(usefcts, fg.IDs[loutn[1].label], fct, fid) #out_neighbors(fct, fg.g)
                    # TODO -- make update vertex call
                    fct.attributes["data"].potentialused = true #fct.attributes["potentialused"] = true
                    dlapi.updatevertex!(fg, fct)
                end
                for sepSearch in loutn # out_neighbors(fct, fg.g)
                    if (fg.IDs[sepSearch.label] == fid)
                        continue # skip the fid itself
                    end
                    sea = findmin(abs(allids-fg.IDs[sepSearch.label]))
                    if sea[1]==0.0
                        appendUseFcts!(usefcts, fg.IDs[sepSearch.label], fct, fid)
                        # usefcts = [usefcts;(fg.IDs[sepSearch.label], fct, fid)]
                        fct.attributes["data"].potentialused = true #fct.attributes["potentialused"] = true
                        dlapi.updatevertex!(fg, fct)
                    end
                end
            end
        end
        cliq.attributes["potentials"]=[cliq.attributes["potentials"];usefcts]
    end
    return nothing
end

function getCliquePotentials!(fg::FactorGraph, bt::BayesTree, chkcliq::Int64)
    getCliquePotentials!(fg, bt.cliques[chkcliq])
end

function cliqPotentialIDs(cliq::Graphs.ExVertex)
  potIDs = Int[]
  for idfct in cliq.attributes["potentials"]
    push!(potIDs,idfct[1])
  end
  return potIDs
end

function collectSeparators(bt::BayesTree, cliq::Graphs.ExVertex)
  allseps = Int[]
  for child in out_neighbors(cliq, bt.bt)#tree
    # if length(child.attributes["conditIDs"]) > 0
      allseps = [allseps; child.attributes["conditIDs"]]
    # end
  end
  return allseps
end

function compCliqAssocMatrices!(bt::BayesTree, cliq::Graphs.ExVertex)
  frtl = cliq.attributes["frontalIDs"]
  cond = cliq.attributes["conditIDs"]
  inmsgIDs = collectSeparators(bt, cliq)
  potIDs = cliqPotentialIDs(cliq)
  # Construct associations matrix here
  # matrix has variables are columns, and messages/constraints as rows
  cols = [frtl;cond]
  cliq.attributes["inmsgIDs"] = inmsgIDs
  cliq.attributes["potIDs"] = potIDs
  cliqAssocMat = Array{Bool,2}(length(potIDs), length(cols))
  cliqMsgMat = Array{Bool,2}(length(inmsgIDs), length(cols))
  fill!(cliqAssocMat, false)
  fill!(cliqMsgMat, false)
  for j in 1:length(cols)
    for i in 1:length(inmsgIDs)
      if cols[j] == inmsgIDs[i]
        cliqMsgMat[i,j] = true
      end
    end
    for i in 1:length(potIDs)
      idfct = cliq.attributes["potentials"][i]
      if idfct[2].index == potIDs[i] # sanity check on clique potentials ordering
        for vertidx in idfct[2].attributes["data"].fncargvID #for vert in idfct[2].attributes["fnc"].Xi
          if vertidx == cols[j]
            cliqAssocMat[i,j] = true
          end
        end
        # for vert in idfct[2].attributes["data"].fnc.Xi #for vert in idfct[2].attributes["fnc"].Xi
        #   if vert.index == cols[j]
        #     cliqAssocMat[i,j] = true
        #   end
        # end
      else
        prtslperr("compCliqAssocMatrices! -- potential ID ordering was lost")
      end
    end
  end
  cliq.attributes["cliqAssocMat"] = cliqAssocMat
  cliq.attributes["cliqMsgMat"] = cliqMsgMat
  nothing
end

function getCliqMat(cliq::Graphs.ExVertex; showmsg=true)
  assocMat = cliq.attributes["cliqAssocMat"]
  msgMat = cliq.attributes["cliqMsgMat"]
  mat = showmsg ? [assocMat;msgMat] : assocMat
  # @show size(mat), size(assocMat), size(msgMat)
  return mat
end

function spyCliqMat(cliq::Graphs.ExVertex; showmsg=true)
  mat = getCliqMat(cliq, showmsg=showmsg)
  Gadfly.spy(mat)
end

function countSkips(bt::BayesTree)
  skps = 0
  for cliq in bt.cliques
    m = getCliqMat(cliq[2])
    mi = map(Int,m)
    skps += sum(map(Int,sum(mi,1) .== 1))
  end
  return skps
end

function skipThroughMsgsIDs(cliq::Graphs.ExVertex)
  numfrtl1 = floor(Int,length(cliq.attributes["frontalIDs"])+1)
  condAssocMat = cliq.attributes["cliqAssocMat"][:,numfrtl1:end]
  condMsgMat = cliq.attributes["cliqMsgMat"][:,numfrtl1:end]
  mat = [condAssocMat;condMsgMat];
  mab = sum(map(Int,mat),1) .== 1
  mabM = sum(map(Int,condMsgMat),1) .== 1
  mab = mab & mabM
  # rang = 1:size(condMsgMat,2)
  msgidx = cliq.attributes["conditIDs"][collect(mab)]
  return msgidx
end

function directAssignmentIDs(cliq::Graphs.ExVertex)
  assocMat = cliq.attributes["cliqAssocMat"]
  msgMat = cliq.attributes["cliqMsgMat"]
  mat = [assocMat;msgMat];
  mab = sum(map(Int,mat),1) .== 1
  mabA = sum(map(Int,assocMat),1) .== 1
  mab = mab & mabA
  frtl = cliq.attributes["frontalIDs"]
  cond = cliq.attributes["conditIDs"]
  cols = [frtl;cond]
  return cols[collect(mab)]
end

function directFrtlMsgIDs(cliq::Graphs.ExVertex)
  numfrtl = length(cliq.attributes["frontalIDs"])
  frntAssocMat = cliq.attributes["cliqAssocMat"][:,1:numfrtl]
  frtlMsgMat = cliq.attributes["cliqMsgMat"][:,1:numfrtl]
  mat = [frntAssocMat; frtlMsgMat];
  mab = sum(map(Int,mat),1) .== 1
  mabM = sum(map(Int,frtlMsgMat),1) .== 1
  mab = mab & mabM
  return cliq.attributes["frontalIDs"][collect(mab)]
end

function mcmcIterationIDs(cliq::Graphs.ExVertex)
  assocMat = cliq.attributes["cliqAssocMat"]
  msgMat = cliq.attributes["cliqMsgMat"]
  mat = [assocMat;msgMat];
  sum(sum(map(Int,mat),1)) == 0 ? error("mcmcIterationIDs -- unaccounted variables") : nothing
  mab = 1 .< sum(map(Int,mat),1)
  frtl = cliq.attributes["frontalIDs"]
  cond = cliq.attributes["conditIDs"]
  cols = [frtl;cond]
  return cols[collect(mab)]
end

function setCliqMCIDs!(cliq::Graphs.ExVertex)
  cliq.attributes["itervarIDs"] = mcmcIterationIDs(cliq)
  cliq.attributes["msgskipIDs"] = skipThroughMsgsIDs(cliq)
  cliq.attributes["directvarIDs"] = directAssignmentIDs(cliq)
  cliq.attributes["directFrtlMsgIDs"] = directFrtlMsgIDs(cliq)
  nothing
end


# function drawCliqAssocMat(fgl::FactorGrapha, cliq::Graphs.ExVertex)
#   cliqAssocMat = cliq.attributes["cliqAssocMat"]
#   ro, co = size(cliqAssocMat)
#   [@show fg.f[i].label for i in tree.cliques[1].attributes["potIDs"]];
#
#   showmat = Array{ASCIIString
#
# end

# post order tree traversal and build potential functions
function buildCliquePotentials(fg::FactorGraph, bt::BayesTree, cliq::Graphs.ExVertex)
    for child in out_neighbors(cliq, bt.bt)#tree
        buildCliquePotentials(fg, bt, child)
    end
    println("Get potentials $(cliq.attributes["label"])");
    getCliquePotentials!(fg, bt, cliq);

    compCliqAssocMatrices!(bt, cliq);
    setCliqMCIDs!(cliq);

    return Union{}
end
