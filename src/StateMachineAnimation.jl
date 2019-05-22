
getStateLabel(state) = Symbol(split(split(string(state), '_')[1],'.')[end])


function histStateMachineTransitions(hist::T) where {T <: Array}
  # local memory
  stateVisits = Dict{Symbol, Vector{Symbol}}()
  allStates = Vector{Symbol}()

  # find all states and transitions
  for i in 1:(length(hist)-1)
    sta  = string(hist[i][3])
    lbl = getStateLabel(sta)
    nsta = string(hist[i+1][3])
    nlbl = getStateLabel(nsta)
    if !haskey(stateVisits, lbl)
      stateVisits[lbl] = Symbol[nlbl;]
    else
      push!(stateVisits[lbl], nlbl)
    end
    union!(allStates, [lbl; nlbl])
  end

  return stateVisits, allStates
end


function histGraphStateMachineTransitions(stateVisits, allStates)

  g = Graphs.incdict(Graphs.ExVertex,is_directed=true)
  lookup = Dict{Symbol, Int}()

  fid = 0
  for state in allStates
    fid += 1
    lbl = string(state)
    exvert = Graphs.ExVertex(fid, lbl)
    exvert.attributes["label"] = lbl
    Graphs.add_vertex!(g, exvert)
    lookup[state] = fid
  end

  count = 0
  for (from, tos) in stateVisits
    for to in tos
      count += 1
      exvf = g.vertices[lookup[from]]
      exvt = g.vertices[lookup[to]]
      # add the edge fom one to the next state
      edge = Graphs.make_edge(g, exvf, exvt)
      Graphs.add_edge!(g, edge)
    end
  end

  return g, lookup
end



function drawStateTransitionStep(hist,
                                 step::Int,
                                 vg,
                                 lookup::Dict{Symbol,Int};
                                 title::String="",
                                 viewerapp="eog",
                                 fext="png",
                                 engine = "dot",
                                 show::Bool=true,
                                 folder::String="",
                                 frame::Int=step  )
  #

  mkpath("/tmp/$folder/")
  dotfile = "/tmp/$folder/csm_$frame.dot"
  filepath="/tmp/$folder/csm_$frame.$(fext)"

  for (id, vert) in vg.vertices
    delete!(vert.attributes, "fillcolor")
    delete!(vert.attributes, "style")
  end

  lbl = getStateLabel(hist[step][3])
  vertid = lookup[lbl]
  vg.vertices[vertid].attributes["fillcolor"] = "red"
  vg.vertices[vertid].attributes["style"] = "filled"

  fid = open(dotfile,"w")
  write(fid,Graphs.to_dot(vg))
  close(fid)

  timest = split(string(hist[step][1]),'T')[end]
  fid = open("/tmp/$folder/dotscript.sh","w")
  str = "head -n `wc -l $dotfile | awk '{print \$1-1}'` $dotfile > /tmp/$folder/tmpdot.dot"
  println(fid, str)
  println(fid, "echo \"graph [label=\\\"$title, #$step, $(timest)\\\", labelloc=t];\" >> /tmp/$folder/tmpdot.dot")
  println(fid, "echo \"}\" >> /tmp/$folder/tmpdot.dot")
  close(fid)
  run(`chmod u+x /tmp/$folder/dotscript.sh`)
  run(`sh /tmp/$folder/dotscript.sh`)
  Base.rm(dotfile)
  Base.rm("/tmp/$folder/dotscript.sh")
  run(`mv /tmp/$folder/tmpdot.dot $dotfile`)

  run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)
  show ? (@async @async run(`$viewerapp $filepath`)) : nothing

  return filepath
end





function drawStateMachineHistory(hist; show::Bool=false, folder::String="")

  stateVisits, allStates = histStateMachineTransitions(hist)

  vg, lookup = histGraphStateMachineTransitions(stateVisits, allStates)

  for i in 1:length(hist)
    drawStateTransitionStep(hist, i, vg, lookup, folder=folder, show=show)
  end

  return nothing
end




function animateStateMachineHistoryByTime(hist::Vector{Tuple{DateTime, Int, <: Function, T}};
                                          frames::Int=100,
                                          folder="animatestate",
                                          title::String="",
                                          show::Bool=false,
                                          startT=hist[1][1],
                                          stopT=hist[end][1]  ) where T
  #
  stateVisits, allStates = histStateMachineTransitions(hist)

  vg, lookup = histGraphStateMachineTransitions(stateVisits, allStates)

  totT = stopT - startT

  step = 1
  len = length(hist)
  @showprogress "exporting state machine images, $title " for i in 1:frames
    aniT = i/frames*totT + startT
    if hist[step][1] < aniT && step < len
      step += 1
    end
    drawStateTransitionStep(hist, step, vg, lookup, title=title, folder=folder, show=show, frame=i)
  end

  nothing
end
