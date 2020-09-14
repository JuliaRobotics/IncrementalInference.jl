
export CSMOccuranceType
export parseCSMVerboseLog, calcCSMOccurancesFolders, calcCSMOccuranceMax, printCSMOccuranceMax, reconstructCSMHistoryLogical

# [cliqId][fsmIterNumber][fsmFunctionName] => (nr. call occurances, list global call sequence position)
const CSMOccuranceType = Dict{Int, Dict{Int, Dict{Symbol, Tuple{Int, Vector{Int}}}}}

function parseCSMVerboseLog(resultsDir::AbstractString;verboseName::AbstractString="csmVerbose.log")
  #
  fid = open(joinpath(resultsDir, verboseName), "r")
  fsmLines = readlines(fid)
  close(fid)

  # parse lines into usable format
  sfsmL = split.(fsmLines, r" -- ")
  cids = split.(sfsmL .|> x->match(r"cliq\d+", x[1]).match, r"cliq") .|> x->parse(Int,x[end])
  iters = split.(sfsmL .|> x->match(r"iter=\d+", x[1]).match, r"iter=") .|> x->parse(Int,x[end])
  smfnc = sfsmL .|> x->x[2] .|> Symbol

  return cids, iters, smfnc
end


## Make lookup from all runs
function calcCSMOccurancesFolders(folderList::Vector{<:AbstractString};
                                  verboseName::AbstractString="csmVerbose.log" )
  #
  # lookup for histogram on each step per fsm
  # [cliqId][fsmIterNumber][fsmFunctionName] => (nr. call occurances, list global call sequence position)
  csmCounter = CSMOccuranceType()
  # lookup for transition counts per fsm function
  trxCounter = Dict{Symbol, Dict{Symbol, Int}}()

  prevFnc = Dict{Int, Symbol}()

  for rDir in folderList
    ## load the sequence from each file
    cids, iters, smfnc = parseCSMVerboseLog(rDir, verboseName=verboseName)

    # populate histogram
    for (idx,smfi) in enumerate(smfnc)
      if !haskey(csmCounter, cids[idx])
        csmCounter[cids[idx]] = Dict{Int, Dict{Symbol, Tuple{Int,Vector{Int}}}}()
      end
      if !haskey(csmCounter[cids[idx]], iters[idx])
        # Tuple{Int,Int[]} == (nr. occurances of call, list global call sequence position)
        csmCounter[cids[idx]][iters[idx]] = Dict{Symbol, Tuple{Int,Vector{Int}}}()
      end
      easyRef = csmCounter[cids[idx]][iters[idx]]
      if !haskey(easyRef,smfi)
        easyRef[smfi] = (0,Int[])
      end
      # add position in call sequence (global per solve)
      globalSeqIdx = easyRef[smfi][2]
      push!(globalSeqIdx, idx)
      easyRef[smfi] = (easyRef[smfi][1]+1, globalSeqIdx)

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



"""
    $SIGNATURES

Use maximum occurance from `csmCounter::CSMOccuranceType` to summarize many CSM results.

Notes
- `percentage::Bool=false` shows median global sequence occurance ('m'), or 
  - `percentage::Bool=true` of occurance ('%')
"""
function calcCSMOccuranceMax( csmCounter::CSMOccuranceType;
                              percentage::Bool=false)
  #
  ncsm = length(keys(csmCounter))
  maxOccuran = Dict()
  # max steps
  for i in 1:ncsm
    # sequence of functions that occur most often
    maxOccuran[i] = Vector{Pair{Symbol, String}}() 
  end

  # pick out the max for each CSM iter
  for (csmID, csmD) in csmCounter, stp in 1:length(keys(csmD))
    maxFnc = :null
    maxCount = 0
    totalCount = 0
    for (fnc, cnt) in csmCounter[csmID][stp]
      totalCount += cnt[1]
      if maxCount < cnt[1]
        maxCount = cnt[1]
        maxFnc = fnc
      end
    end
    perc = if percentage
      "$(round(Int,(maxCount/totalCount)*100))"
    else
      "$(round(Int,Statistics.median(csmCounter[csmID][stp][maxFnc][2])))"
    end
    push!(maxOccuran[csmID], maxFnc=>perc) # position in vector == stp
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
                              fid=stdout,
                              percentage::Bool=false )
  #
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
        fncName = maxOcc[cid][stp][1]
        # either show percentage or sequence index
        percOrSeq = "$(maxOcc[cid][stp][2])"
        percOrSeq *= percentage ? "%" : "m"
        tpl = ("",percOrSeq,fncName,"")
      end
      push!(TPL, tpl)
    end
    IIF.printHistoryLane(fid, stp, TPL)
  end
end



"""
    $SIGNATURES

Use `solveTree!`'s` `verbose` output to reconstruct the swim lanes Logical sequence of CSM function calls.

Notes
- This is a secondary function to primary `printCSMHistoryLogical`.

Related

printCSMHistoryLogical
"""
function reconstructCSMHistoryLogical(resultsDir::AbstractString;
                                      fid::IO=stdout,
                                      verboseName::AbstractString="csmVerbose.log" )
  #
  csmCounter, trxCounter = calcCSMOccurancesFolders([resultsDir], verboseName=verboseName)

  # print with sequence position
  maxOcc = calcCSMOccuranceMax(csmCounter, percentage=false)

  printCSMOccuranceMax(maxOcc, fid=fid)
end






#