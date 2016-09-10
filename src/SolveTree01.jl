

type EasyMessage
  pts::Array{Float64,2}
  bws::Array{Float64,1}
end

type NBPMessage <: Singleton
  p::Dict{Int64,EasyMessage}
end
# type NBPMessage <: Singleton
#     Xi::Array{Int64,1}
#     p::Array{BallTreeDensity,1}
# end

type PotProd
    Xi::Int64
    prev::Array{Float64,2}
    product::Array{Float64,2}
    potentials::Array{BallTreeDensity,1}
    # PotProd() = new()
    # PotProd(x...) = new(x[1],x[2],x[3])
end
type CliqGibbsMC
    prods::Array{PotProd,1}
    # CliqGibbsMC() = new()
    # CliqGibbsMC(x) = new(x)
end
type DebugCliqMCMC
    mcmc::Array{CliqGibbsMC,1}
    outmsg::NBPMessage
    # DebugCliqMCMC() = new()
    # DebugCliqMCMC(x...) = new(x[1],x[2])
end

type UpReturnBPType
    upMsgs::NBPMessage #Dict{Int64, Array{Float64,2}}#
    dbgUp::DebugCliqMCMC
    IDvals::Dict{Int64, Array{Float64,2}}
end

type DownReturnBPType
    dwnMsg::NBPMessage
    dbgDwn::DebugCliqMCMC
    IDvals::Dict{Int64,Array{Float64,2}}
end

type ExploreTreeType
  fg::FactorGraph
  bt::BayesTree
  cliq::Graphs.ExVertex
  prnt::Any
  sendmsgs::Array{NBPMessage,1}
end

type MsgPassType
  fg::FactorGraph
  cliq::Graphs.ExVertex
  vid::Int
  msgs::Array{NBPMessage,1} #dwnMsgs
  N::Int
end

# type UpMsgPassType
#   fg::FactorGraph
#   cliq::Graphs.ExVertex
#   upMsgs::Array{NBPMessage,1}
#   N::Int
# end


#global pidx
pidx = 1
pidl = 1
pidA = 1
thxl = nprocs() > 4 ? floor(Int,nprocs()*0.333333333) : 1

# upploc to control processes done local to this machine and separated from other
# highly loaded processes. upploc() should be used for dispatching short and burst
# of high bottle neck computations. Use upp2() for general multiple dispatch.
function upploc()
    global pidl, thxl
    N = nprocs()
    pidl = (thxl > 1 ? ((pidl < thxl && pidl != 0 && thxl > 2) ? (pidl+1)%(thxl+1) : 2) : (N == 1 ? 1 : 2))
    return pidl
end

# upp2() may refer to processes on a different machine. Do not associate upploc()
# with processes on a separate computer -- this will be become more complicated when
# remote processes desire their own short fast 'local' burst computations. We
# don't want all the data traveling back and forth for shorter jobs
function upp2()
  global pidx, thxl
  N = nprocs()
  pidx = (N > 1 ? ((pidx < N && pidx != 0) ? (pidx+1)%(N+1) : thxl+1) : 1) #2 -- thxl+1
  return pidx
end

function uppA()
  global pidA
  N = nprocs()
  pidA = (N > 1 ? ((pidA < N && pidA != 0) ? (pidA+1)%(N+1) : 2) : 1) #2 -- thxl+1
  return pidA
end


function packFromIncomingDensities!(dens::Array{BallTreeDensity,1}, vertid::Int64, inmsgs::Array{NBPMessage,1})
    for m in inmsgs
        #@show Xi = keys(m.p)
        for idx in keys(m.p) #1:length(Xi)
            if idx == vertid #Xi[idx]
                pdi = m.p[vertid]
                push!(dens, kde!(pdi.pts, pdi.bws) ) #m.p[idx]
            end
            # TODO -- we can inprove speed of search for inner loop
        end
    end
    nothing
end

# add all potentials associated with this clique and vertid to dens
function packFromLocalPotentials!(fgl::FactorGraph, dens::Array{BallTreeDensity,1}, cliq::Graphs.ExVertex, vertid::Int64, N::Int64)
    for idfct in cliq.attributes["potentials"]
        for vertidx in idfct[2].attributes["data"].fncargvID #for vert in idfct[2].attributes["fnc"].Xi
            if vertidx == vertid
              p = findRelatedFromPotential(fgl, idfct[2], vertid, N)
              push!(dens, p)
            end
        end
        # for vert in idfct[2].attributes["data"].fnc.Xi #for vert in idfct[2].attributes["fnc"].Xi
        #     if vert.index == vertid
        #       p = findRelatedFromPotential(fgl, idfct[2], vertid, N)
        #       push!(dens, p)
        #     end
        # end
    end
    nothing
end

function cliqGibbs(fg::FactorGraph, cliq::Graphs.ExVertex, vertid::Int64, inmsgs::Array{NBPMessage,1}, N::Int=200)
    # several optimizations can be performed in this function TODO
    print("$(dlapi.getvertex(fg,vertid).attributes["label"]) ") # "$(fg.v[vertid].attributes["label"]) "
    #consolidate NBPMessages and potentials
    dens = Array{BallTreeDensity,1}()
    packFromIncomingDensities!(dens, vertid, inmsgs)
    packFromLocalPotentials!(fg, dens, cliq, vertid, N)
    # testing
    # end testing
    potprod = PotProd(vertid, getVal(dlapi.getvertex(fg,vertid)), Array{Float64,2}(), dens) # (fg.v[vertid])

    pGM = Array{Float64,2}()
    if length(dens) > 1
        Ndims = dens[1].bt.dims
        dummy = kde!(rand(Ndims,N),[1.0]);
        print("[$(length(dens)) prod. d$(Ndims),N$(N)],")
        pGM, = prodAppxMSGibbsS(dummy, dens, Union{}, Union{}, 8) #10
        #pGM, = remoteProdAppxMSGibbsS(dummy, dens, Union{}, Union{})
        # sum(abs(pGM))<1e-14 ? error("cliqGibbs -- nothing in pGM") : nothing
    elseif length(dens) == 1
        print("[direct]")
        pGM = reshape(dens[1].bt.centers[(dens[1].bt.dims*Npts(dens[1])+1):end], dens[1].bt.dims, Npts(dens[1]))  #N
    else
        pGM = Array{Float64,2}(0,1)
    end
    potprod.product = pGM
    print(", ")

    return pGM, potprod
end

function fmcmc!(fgl::FactorGraph, cliq::Graphs.ExVertex, fmsgs::Array{NBPMessage,1}, IDs::Array{Int64,1}, N::Int64, MCMCIter::Int)
    println("------------- functional mcmc ----------------$(cliq.attributes["label"])")
    # repeat several iterations of functional Gibbs sampling for fixed point convergence
    if length(IDs) == 1
        MCMCIter=1
    end
    mcmcdbg = Array{CliqGibbsMC,1}()

    for iter in 1:MCMCIter
        # iterate through each of the variables, KL-divergence tolerence would be nice test here
        print("#$(iter)\t -- ")
        dbg = CliqGibbsMC([])
        for vertid in IDs
          # we'd like to do this more pre-emptive and then just execute -- just point and skip up only msgs
            densPts, potprod = cliqGibbs(fgl, cliq, vertid, fmsgs, N) #cliqGibbs(fg, cliq, vertid, fmsgs, N)
            if size(densPts,1)>0
                updvert = dlapi.getvertex(fgl,vertid)
                setValKDE!(updvert, densPts) # fgl.v[vertid]
                # Go update the datalayer TODO -- excessive for general case
                dlapi.updatevertex!(fgl, updvert)
                # fgl.v[vertid].attributes["val"] = densPts
                push!(dbg.prods, potprod)
            end
        end
        push!(mcmcdbg, dbg)
        println("")
    end

    # populate dictionary for return NBPMessage in multiple dispatch
    # TODO -- change to EasyMessage dict
    d = Dict{Int64,Array{Float64,2}}()
    for vertid in IDs
        d[vertid] = getVal(dlapi.getvertex(fgl,vertid)) # fgl.v[vertid]
    end
    println("fmcmc! -- finished on $(cliq.attributes["label"])")

    return mcmcdbg, d
end

# function compOutMsgsProd(mpt::MsgPassType)
#     # TODO -- very very expensive function since all fg data is copied on remotecall
#     compIDVals = Array{Float64,2}()
#     compIDVals, = cliqGibbs(mpt.fg, mpt.cliq, mpt.vid, mpt.msgs, mpt.N)
#     if size(compIDVals,1) == 0
#         compIDVals = getVal(mpt.fg.v[mpt.vid])
#     end
#     p = kde!(compIDVals, "lcv")
#     bws = vec(getBW(p)[:,1])
#     return EasyMessage(compIDVals, bws)
# end
#
# # returns output NBPMessage, based on projection to IDs
# function upPrepOutMsg!(fgl::FactorGraph, cliq::Graphs.ExVertex,
#                   inMsgs::Array{NBPMessage,1}, IDs::Array{Int64,1}, N::Int)
#
#     print("Outgoing msg density on: ")
#     len = length(IDs)
#     m = NBPMessage(Dict{Int64,EasyMessage}()) #Int64[], BallTreeDensity[]
#   rr = Array{RemoteRef,1}()
#     i = 0
#     for vid in IDs
#         i+=1
#         mpt = MsgPassType(fgl, cliq, vid, inMsgs, N) # TODO -- should this be copied??
#         if true
#           push!(rr,remotecall(upploc() , compOutMsgsProd, mpt) )
#         else
#           outp = compOutMsgsProd(mpt)
#           m.p[vid] = outp
#         end
#     end
#     i = 0
#     for r in rr
#         i +=1
#         outp = fetch(r)
#         m.p[IDs[i]] = outp
#     end
#     println()
#
#     return m
# end

function upPrepOutMsg!(d::Dict{Int64,Array{Float64,2}}, IDs::Array{Int64,1})
  print("Outgoing msg density on: ")
  len = length(IDs)
  #outDens = Dict{Int64,BallTreeDensity}() # retired to new dictionary msg types
  m = NBPMessage(Dict{Int64,EasyMessage}()) #Int64[], BallTreeDensity[]

  for id in IDs
    p = kde!(d[id], "lcv")
    bws = vec(getBW(p)[:,1])
    m.p[id] = EasyMessage(d[id], bws)
  end
  return m
end


function upGibbsCliqueDensity(inp::ExploreTreeType, N::Int=200)
    print("up w $(length(inp.sendmsgs)) msgs")
    # Loval mcmc over belief functions
    # this is so slow! TODO Can be ignored once we have partial working
    # loclfg = nprocs() < 2 ? deepcopy(inp.fg) : inp.fg

    d = Union{}
    mcmcdbg = Union{}
    # mcmcdbg = [CliqGibbsMC()]

    if false
      IDS = [inp.cliq.attributes["frontalIDs"];inp.cliq.attributes["conditIDs"]] #inp.cliq.attributes["frontalIDs"]
      mcmcdbg, d = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, IDS, N, 3)
    elseif true
      dummy, d = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["directFrtlMsgIDs"], N, 1)
      if length(inp.cliq.attributes["msgskipIDs"]) > 0
        dummy, dd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["msgskipIDs"], N, 1)
        for md in dd d[md[1]] = md[2]; end
      end
      if length(inp.cliq.attributes["itervarIDs"]) > 0
        mcmcdbg, ddd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["itervarIDs"], N, 3)
        for md in ddd d[md[1]] = md[2]; end
      end
      # TODO -- do direct conditionals from msg also, before transits and iterations are done.
      if length(inp.cliq.attributes["directvarIDs"]) > 0
        dummy, dddd = fmcmc!(inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["directvarIDs"], N, 1)
        for md in dddd d[md[1]] = md[2]; end
      end
    end

    #m = upPrepOutMsg!(inp.fg, inp.cliq, inp.sendmsgs, condids, N)
    m = upPrepOutMsg!(d, inp.cliq.attributes["conditIDs"])

    # Copy frontal variables back
    # for id in inp.cliq.attributes["frontalIDs"]
    #     inp.fg.v[id].attributes["val"] = loclfg.v[id].attributes["val"] # inp.
    # end
    # @show getVal(inp.fg.v[1])

    mdbg = DebugCliqMCMC(mcmcdbg, m)
    return UpReturnBPType(m, mdbg, d)
end
# else
#   rr = Array{RemoteRef,1}(4)
#   rr[1] = remotecall(upp2(), fmcmc!, inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["directFrtlMsgIDs"], N, 1)
#   rr[2] = remotecall(upp2(), fmcmc!, inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["msgskipIDs"], N, 1)
#   rr[3] = remotecall(upp2(), fmcmc!, inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["itervarIDs"], N, 3)
#   rr[4] = remotecall(upp2(), fmcmc!, inp.fg, inp.cliq, inp.sendmsgs, inp.cliq.attributes["directvarIDs"], N, 1)
#   dummy, d = fetch(rr[1])
#   dummy, dd = fetch(rr[2])
#   mcmcdbg, ddd = fetch(rr[3])
#   dummy, dddd = fetch(rr[4])
#
#   for md in dd d[md[1]] = md[2]; end
#   for md in ddd d[md[1]] = md[2]; end
#   for md in dddd d[md[1]] = md[2]; end

function dwnPrepOutMsg(fg::FactorGraph, cliq::Graphs.ExVertex, dwnMsgs::Array{NBPMessage,1}, d::Dict{Int64, Array{Float64,2}})
    # pack all downcoming conditionals in a dictionary too.
    if cliq.index != 1
      println("Dwn msg keys $(keys(dwnMsgs[1].p))")
    end # ignore root, now incoming dwn msg
    print("Outgoing msg density on: ")
    # cdwndict = Dict{Int64, BallTreeDensity}()
    # for cvid in cliq.attributes["conditIDs"] # root has no parent
    #     prtdwnmsg = dwnMsgs[1] # there is only one dwn incoming NBPMessage
    #     Xi = keys(prtdwnmsg)
    #     for i in 1:length(Xi) # prtdwnmsg.Xi
    #         if Xi[i] == cvid
    #             print("$(fg.v[cvid].attributes["label"]), ")
    #             p = prtdwnmsg.p[cvid]
    #             cdwndict[cvid] = kde!(p.pts, p.bws)
    #         end
    #     end
    # end
    # println()
    # allIDs = [cliq.attributes["frontalIDs"]; cliq.attributes["conditIDs"] ]
    # outDens = Array{BallTreeDensity,1}(length(allIDs))
    m = NBPMessage(Dict{Int64,EasyMessage}())
    i = 0
    for vid in cliq.attributes["frontalIDs"]
        outp = kde!(d[vid], "lcv") # need to find new down msg bandwidths
        # i+=1
        # outDens[i] = outp
        bws = vec((getBW(outp))[:,1])
        m.p[vid] = EasyMessage(  deepcopy(d[vid]) , bws  )
    end
    for cvid in cliq.attributes["conditIDs"]
        i+=1
        # TODO -- convert to points only since kde replace by rkhs in future
        # outDens[i] = cdwndict[cvid]
        println("")
        println("Looking for cvid=$(cvid)")
        m.p[cvid] = deepcopy(dwnMsgs[1].p[cvid]) # TODO -- maybe this can just be a union(,)
    end

    # NBPMessage(deepcopy(allIDs), deepcopy(outDens))
    return m
end

function downGibbsCliqueDensity(fg::FactorGraph, cliq::Graphs.ExVertex, dwnMsgs::Array{NBPMessage,1}, N::Int=200, MCMCIter::Int=3)
    print("dwn")
    mcmcdbg, d = fmcmc!(fg, cliq, dwnMsgs, cliq.attributes["frontalIDs"], N, MCMCIter)
    m = dwnPrepOutMsg(fg, cliq, dwnMsgs, d)

    mdbg = DebugCliqMCMC(mcmcdbg, m)
    return DownReturnBPType(m, mdbg, d)
end


function updateFGBT!(fg::FactorGraph, bt::BayesTree, cliqID::Int64, ddt::DownReturnBPType)
    # if dlapi.cgEnabled
    #   return nothing
    # end
    cliq = bt.cliques[cliqID]
    # cliq.attributes["debugDwn"] = deepcopy(ddt.dbgDwn) #inp.
    for dat in ddt.IDvals
      #TODO -- should become an update call
        updvert = dlapi.getvertex(fg,dat[1])
        setValKDE!(updvert, deepcopy(dat[2])) # TODO -- not sure if deepcopy is required
        updvert.attributes["latestEst"] = Base.mean(dat[2],2)
        # fg.v[dat[1]].attributes["val"] = deepcopy(dat[2]) # inp.
        dlapi.updatevertex!(fg, updvert)
    end
    nothing
end

# TODO -- use Union{} for two types, rather than separate functions
function updateFGBT!(fg::FactorGraph, bt::BayesTree, cliqID::Int64, urt::UpReturnBPType)
    # if dlapi.cgEnabled
    #   return nothing
    # end
    cliq = bt.cliques[cliqID]
    # cliq.attributes["debug"] = deepcopy(urt.dbgUp) #inp.
    for dat in urt.IDvals
      updvert = dlapi.getvertex(fg,dat[1])
      setValKDE!(updvert, deepcopy(dat[2])) # (fg.v[dat[1]], ## TODO -- not sure if deepcopy is required
      updvert.attributes["latestEst"] = Base.mean(dat[2],2)
      # fg.v[dat[1]].attributes["val"] = deepcopy(dat[2]) # inp.
      dlapi.updatevertex!(fg, updvert)
    end
    println("updateFGBT! up -- finished updating $(cliq.attributes["label"])")
    nothing
end

# pass NBPMessages back down the tree -- pre order tree traversal
function downMsgPassingRecursive(inp::ExploreTreeType; N::Int=200)
    println("====================== Clique $(inp.cliq.attributes["label"]) =============================")

    mcmciter = inp.prnt != Union{} ? 3 : 0; # skip mcmc in root on dwn pass
    rDDT = downGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, mcmciter) #dwnMsg
    updateFGBT!(inp.fg, inp.bt, inp.cliq.index, rDDT)

    # rr = Array{RemoteRef,1}()
    pcs = procs()

    ddt=Union{}
    for child in out_neighbors(inp.cliq, inp.bt.bt)
        ett = ExploreTreeType(inp.fg, inp.bt, child, inp.cliq, [rDDT.dwnMsg])#inp.fg
        ddt = downMsgPassingRecursive( ett , N=N )
        # updateFG!(inp, ddt)
        #r = remotecall(upp2() , downMsgPassingRecursive, ett) #pidx
        #push!(rr, r)
        ##push!(ETT, ett)
        ##r = @spawn downMsgCliquePotentials( ett ) # loss of type information
    end
    # for r in rr
    #     rDDT = fetch(r)
    #     # calling updateFG here is obsolete and probably wrong at this point
    #     println("dwnMsgPR -- calling updateFGBT! on $(inp.cliq.attributes["label"])...")
    #     updateFGBT!(inp.fg, inp.bt, inp.cliq.index, rDDT) # this runs in the main thread (NBPMessage being passed back to proc1)
    # end

    # return modifications to factorgraph to calling process
    # TODO -- is this the right information to be passing down?
    return ddt
end



# post order tree traversal and build potential functions
function upMsgPassingRecursive(inp::ExploreTreeType; N::Int=200) #upmsgdict = Dict{Int64, Array{Float64,2}}()
    println("Start Clique $(inp.cliq.attributes["label"]) =============================")
    childMsgs = Array{NBPMessage,1}()

    len = length(out_neighbors(inp.cliq, inp.bt.bt))
    for child in out_neighbors(inp.cliq, inp.bt.bt)
        ett = ExploreTreeType(inp.fg, inp.bt, child, inp.cliq, NBPMessage[]) # ,Union{})
        println("upMsgRec -- calling new recursive on $(ett.cliq.attributes["label"])")
        # newmsgs = Dict{Int64, Array{Float64,2}}()
        newmsgs = upMsgPassingRecursive(  ett, N=N ) # newmsgs
        println("upMsgRec -- finished with $(ett.cliq.attributes["label"]), w $(keys(newmsgs.p)))")
          # Xi = Array{Int64,1}()
          # p = Array{BallTreeDensity,1}()
          # for umsg in newmsgs
          #     push!(Xi, umsg[1])
          #     push!(p, kde!(umsg[2], "lcv"))
          # end
          # push!(  childMsgs, NBPMessage(Xi, p))
        push!(  childMsgs, newmsgs )
    end

    println("====================== Clique $(inp.cliq.attributes["label"]) =============================")
    ett = ExploreTreeType(inp.fg, inp.bt, inp.cliq, Union{}, childMsgs)

    # @show "BEFORE", getVal(ett.fg.v[1])
    # upmsg = Dict{Int64, Array{Float64,2}}()
    urt = upGibbsCliqueDensity(ett, N) # upmsgdict
    updateFGBT!(inp.fg, inp.bt, inp.cliq.index, urt)
    # @show "AFTER", getVal(ett.fg.v[1])
    println("End Clique $(inp.cliq.attributes["label"]) =============================")
    urt.upMsgs
end


function downGibbsCliqueDensity(inp::ExploreTreeType, N::Int=200)
  println("=================== Iter Clique $(inp.cliq.attributes["label"]) ===========================")
  mcmciter = inp.prnt != Union{} ? 3 : 0
  return downGibbsCliqueDensity(inp.fg, inp.cliq, inp.sendmsgs, N, mcmciter)
end

function prepDwnPreOrderStack!(bt::BayesTree, parentStack::Array{Graphs.ExVertex,1})
  # dwn message passing function
  nodedata = Union{}
  tempStack = Array{Graphs.ExVertex,1}()
  push!(tempStack, parentStack[1])

  while ( length(tempStack) != 0 || nodedata != Union{} )
    if nodedata != Union{}
      for child in out_neighbors(nodedata, bt.bt) #nodedata.cliq, nodedata.bt.bt
          ## TODO -- don't yet have the NBPMessage results to propagate belief
          ## should only construct ett during processing
          push!(parentStack, child) #ett
          push!(tempStack, child) #ett
      end
      nodedata = Union{} # its over but we still need to complete computation in each of the leaves of the tree
    else
      nodedata = tempStack[1]
      deleteat!(tempStack, 1)
    end
  end
  nothing
end

function findVertsAssocCliq(fgl::FactorGraph, cliq::Graphs.ExVertex)

  IDS = [cliq.attributes["frontalIDs"];cliq.attributes["conditIDs"]] #inp.cliq.attributes["frontalIDs"]


  error("findVertsAssocCliq -- not completed yet")
  nothing
end

function partialExploreTreeType(pfg::FactorGraph, pbt::BayesTree, cliqCursor::Graphs.ExVertex, prnt, pmsgs::Array{NBPMessage,1})
    # println("starting pett")
    # TODO -- expand this to grab only partial subsection from the fg and bt data structures



    if length(pmsgs) < 1
      return ExploreTreeType(pfg, pbt, cliqCursor, prnt, NBPMessage[])
    else
      return ExploreTreeType(pfg, pbt, cliqCursor, prnt, pmsgs)
    end
    nothing
end

function dispatchNewDwnProc!(fg::FactorGraph, bt::BayesTree, parentStack::Array{Graphs.ExVertex,1}, stkcnt::Int64, refdict::Dict{Int64,RemoteRef}; N::Int=200)
  #ddt = visitNode(nodedata[2])
  #updateFG!(nodedata[2], ddt)
  cliq = parentStack[stkcnt]
  # while !haskey(cliq.attributes, "dwnremoteref") # nodedata.cliq
  while !haskey(refdict, cliq.index) # nodedata.cliq
    #println("Sleeping on lack of remoteref")
    sleep(0.25)
  end

  # rDDT = fetch(cliq.attributes["dwnremoteref"]) #nodedata.cliq
  # delete!(cliq.attributes, "dwnremoteref") # nodedata
  rDDT = fetch(refdict[cliq.index]) #nodedata.cliq
  delete!(refdict,cliq.index) # nodedata

  if rDDT != Union{}
    updateFGBT!(fg, bt, cliq.index, rDDT)
  end

  emptr = BayesTree(Union{}, 0, Dict{Int,Graphs.ExVertex}(), Dict{ASCIIString,Int64}());

  for child in out_neighbors(cliq, bt.bt) # nodedata.cliq, nodedata.bt.bt
      #haskey(child.attributes, "dwnremoteref") ? error("dispatchNewDwnProc! -- why you already have dwnremoteref?") : nothing
      haskey(refdict, child.index) ? error("dispatchNewDwnProc! -- why you already have dwnremoteref?") : nothing

      ett = partialExploreTreeType(fg, emptr, child, cliq, [rDDT.dwnMsg]) # bt
      #child.attributes["dwnremoteref"] = remotecall(upp2() , downGibbsCliqueDensity, ett, N) # pidxI[1]
      refdict[child.index] = remotecall(upp2() , downGibbsCliqueDensity, ett, N) # pidxI[1]
  end
  nothing
end

function processPreOrderStack!(fg::FactorGraph, bt::BayesTree, parentStack::Array{Graphs.ExVertex,1}, refdict::Dict{Int64,RemoteRef}; N::Int=200)
    # dwn message passing function for iterative tree exploration
    # @show length(parentStack)
    stkcnt = 0

    @sync begin
      while ( stkcnt < length(parentStack) ) # || nodedata != Union{}
          stkcnt += 1
          @async dispatchNewDwnProc!(fg, bt, parentStack, stkcnt, refdict, N=N) #pidxI,nodedata
      end
    end
    nothing
end

function downMsgPassingIterative!(startett::ExploreTreeType; N::Int=200)
  # this is where we launch the downward iteration process from
  parentStack = Array{Graphs.ExVertex,1}()
  refdict = Dict{Int64,RemoteRef}()

  # start at the given clique in the tree -- shouldn't have to be the root.
  # haskey(startett.cliq.attributes, "dwnremoteref") ? error("downMsgPassingIterative! -- why you already have dwnremoteref?") : nothing
  pett = partialExploreTreeType(startett.fg, startett.bt, startett.cliq,
                                        startett.prnt, startett.sendmsgs)
  # startett.cliq.attributes["dwnremoteref"] = remotecall(upp2(), downGibbsCliqueDensity, pett, N) #startett)
  refdict[startett.cliq.index] = remotecall(upp2(), downGibbsCliqueDensity, pett, N)

  push!(parentStack, startett.cliq ) # r

  prepDwnPreOrderStack!(startett.bt, parentStack)
  processPreOrderStack!(startett.fg, startett.bt, parentStack, refdict, N=N)

  println("dwnward leftovers, $(keys(refdict))")
  # println("dwnStk counter at end is $(stkcnt), length(parentStack) = $(length(parentStack))")
  nothing
end

function prepPostOrderUpPassStacks!(bt::BayesTree, parentStack::Array{Graphs.ExVertex,1}, childStack::Array{Graphs.ExVertex,1})
  # upward message passing preparation
  while ( length(parentStack) != 0 )
      #2.1 Pop a node from first stack and push it to second stack
      cliq = parentStack[end]
      deleteat!(parentStack, length(parentStack))
      push!(childStack, cliq )

      #2.2 Push left and right children of the popped node to first stack
      for child in out_neighbors(cliq, bt.bt) # nodedata.cliq, nodedata.bt.bt
          push!(parentStack, child )
      end
  end
  for child in childStack
    @show child.attributes["label"]
  end
  nothing
end

function asyncProcessPostStacks!(fgl::FactorGraph, bt::BayesTree, chldstk::Array{Graphs.ExVertex,1}, stkcnt::Int, refdict::Dict{Int64,RemoteRef};
                                  N::Int=200)
  # for up message passing
  cliq = chldstk[stkcnt]
  gomulti = true
  println("Start Clique $(cliq.attributes["label"]) =============================")
  childMsgs = Array{NBPMessage,1}()
  # println("asyncProcessPostStacks -- new async at stkcnt=$(stkcnt), cliq=$(cliq.attributes["label"]) has $(length(out_neighbors(cliq, bt.bt))) children")
  ur = Union{}
  for child in out_neighbors(cliq, bt.bt) #nodedata.cliq, nodedata.bt.bt
      println("asyncProcessPostStacks -- $(stkcnt), cliq=$(cliq.attributes["label"]), start on child $(child.attributes["label"]) haskey=$(haskey(child.attributes, "remoteref"))")

        #while !haskey(child.attributes, "remoteref")
        while !haskey(refdict, child.index)
          # println("Sleeping $(cliq.attributes["label"]) on lack of remoteref from $(child.attributes["label"])")
          sleep(0.25)
        end

      if gomulti
        #ur = fetch(child.attributes["remoteref"])
        ur = fetch(refdict[child.index])
      else
        ur = child.attributes["remoteref"]
      end
      updateFGBT!( fgl, bt, child.index, ur ) # deep copies happen in the update function
      #delete!(child.attributes, "remoteref")

      push!(childMsgs, ur.upMsgs)

  end
  println("====================== Clique $(cliq.attributes["label"]) =============================")
  #sleep(2.0) # delay sometimes makes it work -- we might have a race condition
  emptr = BayesTree(Union{}, 0, Dict{Int,Graphs.ExVertex}(), Dict{ASCIIString,Int64}());
  pett = partialExploreTreeType(fgl, emptr, cliq, Union{}, childMsgs) # bt   # parent cliq pointer is not needed here, fix Graphs.jl first

  if haskey(cliq.attributes, "remoteref")
      println("asyncProcessPostStacks! -- WHY YOU ALREADY HAVE REMOTEREF?")
  end

  newprocid = upp2()
  #println("asyncProcessPostStacks -- making remote call to $(newprocid) for $(pett.cliq.attributes["label"]) with $(length(pett.sendmsgs))")
  if gomulti
    #cliq.attributes["remoteref"] = remotecall(newprocid, upGibbsCliqueDensity, pett, N )
    refdict[cliq.index] = remotecall(newprocid, upGibbsCliqueDensity, pett, N )
  else
    cliq.attributes["remoteref"] = upGibbsCliqueDensity(pett, N)
  end

  # delete as late as possible, but could happen sooner
  for child in out_neighbors(cliq, bt.bt)
      # delete!(child.attributes, "remoteref")
      delete!(refdict, child.index)
  end

  println("End Clique $(cliq.attributes["label"]) =============================")
  nothing
end

      # Xi = Array{Int64,1}()
      # p = Array{BallTreeDensity,1}()
      # for umsg in ur.upMsgs
      #     push!(Xi, umsg[1])
      #     push!(p, kde!(umsg[2], "lcv"))
      # end
      # push!(  childMsgs, NBPMessage(Xi, p))


function processPostOrderStacks!(fg::FactorGraph, bt::BayesTree, childStack::Array{Graphs.ExVertex,1}, N::Int=200)
  # upward message passing function

  refdict = Dict{Int64,RemoteRef}()

  stkcnt = length(childStack)
  @sync begin
    while (stkcnt > 0) #length(childStack)
        #nodedata = childStack[end]
        #deleteat!(childStack, length(childStack))
        @async asyncProcessPostStacks!(fg, bt, childStack, stkcnt, refdict, N=N)
        stkcnt -= 1
    end
  end
  println("processPostOrderStacks! -- THIS ONLY HAPPENS AFTER SYNC")
  # we still need to fetch the root node computational output
  # println("processPostOrderStacks! -- going to fetch remoteref")
  # @show childStack[1].attributes["remoteref"]
  if true
    # ur = fetch(childStack[1].attributes["remoteref"])
    ur = fetch(refdict[childStack[1].index])
  else
    ur = childStack[1].attributes["remoteref"]
  end
  # delete!(childStack[1].attributes, "remoteref") # childStack[1].cliq
  delete!(refdict, childStack[1].index)

  println("upward leftovers, $(keys(refdict))")

  updateFGBT!(fg, bt, childStack[1].index, ur ) # nodedata
  nothing
end

function upMsgPassingIterative!(startett::ExploreTreeType; N::Int=200)
  #http://www.geeksforgeeks.org/iterative-postorder-traversal/
  # this is where we launch the downward iteration process from
  parentStack = Array{Graphs.ExVertex,1}()
  childStack = Array{Graphs.ExVertex,1}()
  #Loop while first stack is not empty
  push!(parentStack, startett.cliq )
  # Starting at the root means we have a top down view of the tree
  prepPostOrderUpPassStacks!(startett.bt, parentStack, childStack)
  # println("upMsgPassingIterative! -- okay lets work through the childStack")
  processPostOrderStacks!(startett.fg, startett.bt, childStack, N)
  nothing
end


function inferOverTree!(fgl::FactorGraph, bt::BayesTree; N::Int=200)
    println("Do inference over tree")
    # Profile.clear()
    cliq = bt.cliques[1]
    upMsgPassingIterative!(ExploreTreeType(fgl, bt, cliq, Union{}, NBPMessage[]),N=N);
    cliq = bt.cliques[1]
    downMsgPassingIterative!(ExploreTreeType(fgl, bt, cliq, Union{}, NBPMessage[]),N=N);
    nothing
end

function inferOverTreeR!(fgl::FactorGraph, bt::BayesTree; N::Int=200)
    println("Do inference over tree")
    # Profile.clear()
    cliq = bt.cliques[1]
    upMsgPassingRecursive(ExploreTreeType(fgl, bt, cliq, Union{}, NBPMessage[]), N=N);
    cliq = bt.cliques[1]
    downMsgPassingRecursive(ExploreTreeType(fgl, bt, cliq, Union{}, NBPMessage[]), N=N);
    nothing
end


### was downGibbsCliqueDensity
# if length(IDs) == 1
#     MCMCIter=1
# end
#
# mcmcdbg = Array{CliqGibbsMC,1}()
# for iter in 1:MCMCIter
#     dbg = CliqGibbsMC([])
#     # iterate through each of the frontal variables in the clique, KL-divergence tolerence would be better
#     print("#$(iter)\t -- ")
#     for vertid in IDs
#         densPts, potprod = cliqGibbs(fg, cliq, vertid, dwnMsgs, N)
#         if size(densPts,1)>0
#             fg.v[vertid].attributes["val"] = densPts
#             push!(dbg.prods, potprod)
#         end
#     end
#     push!(mcmcdbg, dbg)
#     println("")
# end

# populate dictionary for return NBPMessage in multiple dispatch
# d = Dict{Int64,Array{Float64,2}}()
# for vertid in IDs
#     d[vertid] = fg.v[vertid].attributes["val"]
# end
# print("Outgoing msg density on: ")
# outDens = Array{BallTreeDensity,1}(length(allIDs))
# rr = Array{RemoteRef,1}()
# # this can be in parallel
# for vid in allIDs #IDs
#     #pidx = upp(pidx, nprocs())
#     mpt = MsgPassType(fg, cliq, vid, dwnMsgs, N)
#     if false
#       push!(rr,remotecall(upp2() , compOutMsgsProd, mpt) )
#     else
#       outDens[i] = compOutMsgsProd(mpt)
#     end
# end
# count=0
# for r in rr
#   count +=1
#   outDens[count] = fetch(r)
# end
#
# println()
# #@show "outMsg",cliq.attributes["conditIDs"], size(outDens)
# m=NBPMessage(deepcopy(allIDs), deepcopy(outDens))
