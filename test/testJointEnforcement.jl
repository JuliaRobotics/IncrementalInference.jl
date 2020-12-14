# test case to enforce consistency in joint gibbs


using IncrementalInference


##

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{1})
addVariable!(fg, :x1, ContinuousEuclid{1})
addVariable!(fg, :x2, ContinuousEuclid{1})

initManual!(fg, :x0, randn(1,100))
initManual!(fg, :x1, randn(1,100) .+ 10)
initManual!(fg, :x2, randn(1,100) .+ 20)

addFactor!(fg , [:x0; :x1], LinearRelative(MvNormal([10.0;], diagm([1.0;]))))
addFactor!(fg , [:x1; :x2], LinearRelative(MvNormal([10.0;], diagm([1.0;]))))

# setPPE!(fg, :x2)
# fg[:x2]


##

addVariable!(fg, :x3, ContinuousScalar)
addFactor!(fg, [:x2; :x3], LinearRelative(MvNormal([10.0], diagm([1.0]))))

##

addFactor!(fg, [:x0; :x3], LinearRelative(MvNormal([10.0], diagm([1.0]))), graphinit=false)


##

drawGraph(fg, show=true)

##

tree, _, = solveTree!(fg);


##

drawTree(tree, show=true)

## get up message from child clique

msgBuf = IIF.getMessageBuffer(getClique(tree, :x3))

msg = msgBuf.upTx



## findSubgraphs

  ## new version needs to search better which diffJoints are included
  # 1. count separtor connectivity in UPWARD_DIFFERENTIAL
  # 2. start with 0's as subgraphs
  # 3a. then 1's and search all paths, adding each hit to subgraph lists
  # 3b. then 2's
  # 4. until all separators allocated to subgraph



commonJoints = []
subClassify = Dict{Symbol,Int}()
newClass = 0

separators = [:x2;:x0; :x4]

# 1. count separtor connectivity in UPWARD_DIFFERENTIAL
sepsCount = Dict{Symbol, Int}()
map(x->(sepsCount[x]=0), separators)
# tagsFilter = [:LIKELIHOODMESSAGE;]
# tflsf = lsf(fg, tags=tagsFilter)
for likl in msg.diffJoints
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
for key in setdiff(keys(sepsCount), keys(subClassify))
  if !(key in keys(subClassify))
    newClass += 1
    subClassify[key] = newClass
  end
  # if sepsCount[key] == 1
    # search connectivity throughout remaining variables, some duplicate computation occurring
    for key2 in setdiff(keys(sepsCount), keys(subClassify))
      pth = findShortestPathDijkstra(subfg, key, key2)
      # check if connected to existing subClass
      if 0 == length(pth)
        # not connected, so need new class
        newClass += 1
        subClassify[key2] = newClass
      else
        # is connected, so add existing class of key
        subClassify[key2] = subClassify[key]
      end
    end
  # end
end

# 4. inverse classification dictionary
allClasses = Dict{Int, Vector{Symbol}}()
for (key, cls) in subClassify
  !haskey(allClasses, cls) ? (allClasses[cls] = Symbol[key;]) : union!(allClasses[cls],[key;])
end
@show allClasses


# 5. find best variable of each of allClasses to place MsgPrior

bestCandidate = IIF._calcCandidatePriorBest(subfg, msg, allClasses[1])



  # # add a prior to this variable
  # upcm = IIF.generateMsgPrior(TreeBelief(getVariable(subfg, id)), msg.msgType)
  # push!(commonJoints, upcm)


##

mb = IIF.getMessageBuffer(getClique(tree, :x0))

mb.upRx[2].diffJoints[1][1]
mb.upRx[2].diffJoints[1][2]


##

@show findShortestPathDijkstra(fg, :x0,:x2)

@show findShortestPathDijkstra(fg, :x1,:x3)



##

isPathFactorsHomogeneous(fg, :x0, :x2)


##

fg[:x4]
fg[:x2]

##


# vo = getEliminationOrder(fg)

tree = resetBuildTree!(fg)
# tree = resetBuildTreeFromOrder!(fg, [:x0;:x1;:x3;:x2])

drawTree(tree, show=true)



##


# each clique has a subgraph
cfg2 = buildCliqSubgraph(fg,getClique(tree,:x3))
drawGraph(cfg2, show=true)


## check which path between separators has homogeneous factors


isHom, ftyps = isPathFactorsHomogeneous(fg, :x0, :x2)


_sft = selectFactorType(fg, :x0, :x2) 
sft = _sft()

typeof(sft).name == ftyps[1]

getindex(Main, ftyps[1])


##  dev


##


retlist = addLikelihoodsDifferentialCHILD!([:x2; :x0], cfg2)


retlist[1][2]




##

drawGraph(tfg, show=true)


##

getManifolds(LinearRelative)





#
