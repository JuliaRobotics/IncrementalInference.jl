
export calcCSMOccurancesFolders, calcCSMOccuranceMax, printCSMOccuranceMax



## Make lookup from all runs
function calcCSMOccurancesFolders(folderList,
                                  verboseName::AbstractString="csmVerbose.log" )
  #
  # lookup for histogram on each step per fsm
  csmCounter = Dict()
  # lookup for transition counts per fsm function
  trxCounter = Dict{Symbol, Dict{Symbol, Int}}()

  prevFnc = Dict{Int, Symbol}()

  for rDir in folderList
    ## load the sequence from each file
    fid = open(joinpath(rDir, verboseName), "r")
    fsmLines = readlines(fid)
    close(fid)

    # parse lines into usable format
    sfsmL = split.(fsmLines, r" -- ")
    cids = split.(sfsmL .|> x->match(r"cliq\d+", x[1]).match, r"cliq") .|> x->parse(Int,x[end])
    iters = split.(sfsmL .|> x->match(r"iter=\d+", x[1]).match, r"iter=") .|> x->parse(Int,x[end])
    smfnc = sfsmL .|> x->x[2] .|> Symbol

    # populate histogram
    for (idx,smfi) in enumerate(smfnc)
      if !haskey(csmCounter, cids[idx])
        csmCounter[cids[idx]] = Dict{Int, Dict{Symbol, Int}}()
      end
      if !haskey(csmCounter[cids[idx]], iters[idx])
        csmCounter[cids[idx]][iters[idx]] = Dict{Symbol, Int}()
      end
      if !haskey(csmCounter[cids[idx]][iters[idx]],smfi)
        csmCounter[cids[idx]][iters[idx]][smfi] = 0
      end
      csmCounter[cids[idx]][iters[idx]][smfi] += 1

      ## also track the transitions
      if haskey(prevFnc, cids[idx])
        if !haskey(trxCounter, prevFnc[cids[idx]])
          # add function lookup if not previously seen 
          trxCounter[prevFnc[cids[idx]]] = Dict{Symbol, Int}()
        end
        if !haskey(trxCounter[prevFnc[cids[idx]]], smfi)
          # add previously unseen transition
          trxCounter[prevFnc[cids[idx]]][smfi] = 0
        end
        # from previous to next function
        trxCounter[prevFnc[cids[idx]]][smfi] += 1
      end
      # always update prevFnc register
      prevFnc[cids[idx]] = smfi
    end
  end

  return csmCounter, trxCounter
end




function calcCSMOccuranceMax(csmCounter::Dict)
  ncsm = length(keys(csmCounter))
  maxOccuran = Dict()
  # max steps
  for i in 1:ncsm
    # sequence of functions that occur most often
    maxOccuran[i] = Vector{Pair{Symbol, Float64}}() 
  end

  # pick out the max for each CSM iter
  for (csmID, csmD) in csmCounter, stp in 1:length(keys(csmD))
    maxFnc = :null
    maxCount = 0
    totalCount = 0
    for (fnc, cnt) in csmCounter[csmID][stp]
      totalCount += cnt
      if maxCount < cnt
        maxCount = cnt
        maxFnc = fnc
      end
    end
    push!(maxOccuran[csmID], maxFnc=>(maxCount/totalCount)) # position in vector == stp
  end
  maxOccuran
end


"""
    $SIGNATURES

Print the most likely FSM function at each step per state machine, as swim lanes.

Example

```julia
csmCo = calcCSMOccurancesFolders(resultFolder[maskTrue])
maxOcc = calcCSMOccuranceMax(csmCo)
printCSMOccuranceMax(maxOcc)
```
"""
function printCSMOccuranceMax(maxOcc;
                              fid=stdout)
  ncsm = length(keys(maxOcc))

  # print titles
  titles = Tuple[]
  for cid in 1:ncsm
    tpl = ("","","$cid                   ","")
    push!(titles, tpl)
  end
  IIF.printHistoryLane(fid, "", titles)
  print(fid,"----")
  for i in 1:ncsm
    print(fid,"+----------------")
  end
  println(fid,"")


  maxsteps=0
  for i in 1:ncsm
    maxsteps = maxsteps < length(maxOcc[i]) ? length(maxOcc[i]) : maxsteps
  end

  for stp in 1:maxsteps
    TPL = Tuple[]
    for cid in 1:ncsm
      tpl = ("","","                     ","")
      if stp <= length(maxOcc[cid])
        tpl = ("","$(round(Int,maxOcc[cid][stp][2]*100))%",maxOcc[cid][stp][1],"")
      end
      push!(TPL, tpl)
    end
    IIF.printHistoryLane(fid, stp, TPL)
  end
end





#