# test case to enforce consistency in joint gibbs


using IncrementalInference


##

fg = initfg()

addVariable!(fg, :x0, ContinuousEuclid{2})
addVariable!(fg, :x1, ContinuousEuclid{2})
addVariable!(fg, :x2, ContinuousEuclid{2})

initManual!(fg, :x0, randn(2,100))
initManual!(fg, :x1, randn(2,100) .+ 10)
initManual!(fg, :x2, randn(2,100) .+ 20)

addFactor!(fg , [:x0; :x1], LinearRelative(MvNormal([10.0;10], diagm([1.0;1]))))
addFactor!(fg , [:x1; :x2], LinearRelative(MvNormal([10.0;10], diagm([1.0;1]))))

# setPPE!(fg, :x2)
# fg[:x2]


##

addVariable!(fg, :x3, ContinuousScalar)
addFactor!(fg, [:x2; :x3], LinearRelative(MvNormal([10.0;10], diagm([1.0;1]))))

##

addFactor!(fg, [:x0; :x3], LinearRelative(MvNormal([10.0;10], diagm([1.0;1]))), graphinit=false)


##

drawGraph(fg, show=true)

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

# check which path between separators has homogeneous factors
isHom, ftyps = isPathFactorsHomogeneous(fg, :x0, :x2)


_sft = selectFactorType(fg, :x0, :x2) 
sft = _sft()

typeof(sft).name == ftyps[1]

getindex(Main, ftyps[1])


##  dev


# function addLikelihoodsDifferentialCHILD!(seps::Vector{Symbol}, 
#                                           cliqSubFG::AbstractDFG,
#                                           tfg::AbstractDFG=initfg() )
#   #
#   # return this list
#   retlist = Vector{Tuple{Vector{Symbol},DFG.AbstractRelative}}()


#   # create new local dfg and add all the variables with data
#   for label in seps
#     if !exists(tfg, label)
#       addVariable!(tfg, label, getVariableType(fg, label))
#       @debug "New variable added to subfg" _group=:check_addLHDiff #TODO JT remove debug. 
#     end
#     initManual!(tfg, label, getBelief(fg, label))
#   end

#   # list all variables in order of dimension size
#   alreadylist = Symbol[]
#   listDims = getDimension.(getVariable.(tfg, seps))
#   per = sortperm(listDims, rev=true)
#   listVarDec = seps[per]
#   listVarAcc = reverse(listVarDec)
#   # add all differential factors (without deconvolution values)
#   for sym1_ in listVarDec
#     push!(alreadylist, sym1_)
#     for sym2_ in setdiff(listVarAcc, alreadylist)
#       isHom, ftyps = isPathFactorsHomogeneous(fg, :x0, :x2)
#       # chain of user factors are of the same type
#       if isHom
#         @show _sft = selectFactorType(tfg, sym1_, sym2_) 
#         sft = _sft()
#         # only take factors that are homogeneous with the generic relative
#         if typeof(sft).name == ftyps[1]
#           # assume default helper function # buildFactorDefault(nfactype)
#           afc = addFactor!(tfg, [sym1_;sym2_], sft, graphinit=false, tags=[:DUMMY;])
#           # calculate the general deconvolution between variables
#           pts, = approxDeconv(tfg, afc.label)  # solveFactorMeasurements
#           newBel = manikde!(pts, _sft) # getManifolds(sft)
#           # replace dummy factor with real deconv factor using manikde approx belief measurement
#           fullFct = _sft(newBel)
#           deleteFactor!(tfg, afc.label)
#           push!(retlist, ([sym1_;sym2_], fullFct))
#         end
#       end
#     end
#   end

#   return retlist
# end


##


retlist = addLikelihoodsDifferentialCHILD!([:x2; :x0], cfg2)


retlist[1][2]




##

drawGraph(tfg, show=true)


##

getManifolds(LinearRelative)





#
