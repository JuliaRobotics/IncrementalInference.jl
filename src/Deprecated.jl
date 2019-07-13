
"""
    $SIGNATURES

Draw and show the factor graph `<:AbstractDFG` via system graphviz and pdf app.

Notes
- Should not be calling outside programs.
- Need long term solution
- DFG's `toDotFile` a better solution -- view with `xdot` application.
- also try `engine={"sfdp","fdp","dot","twopi","circo","neato"}`

Future:
- Might be kept with different call strategy since this function so VERY useful!
- Major issue that this function calls an external program such as "evince", which should be
   under user control only.
- Maybe solution is
- ```toDot(fg,file=...); @async run(`xdot file.dot`)```, or
  - ```toDot(fg,file=...); exportPdf(...); @async run(`evince ...pdf`)```.
"""
function writeGraphPdf(fgl::G;
                       viewerapp::String="evince",
                       filepath::AS="/tmp/fg.pdf",
                       engine::AS="neato", #sfdp
                       show::Bool=true ) where {G <: AbstractDFG, AS <: AbstractString}
  #
  @warn "writeGraphPdf function might changed, see DFG.toDotFile(dfg) as part of the long term solution."

  fgd = fgl
  @info "Writing factor graph file"
  fext = split(filepath, '.')[end]
  fpwoext = split(filepath, '.')[end-1]
  dotfile = fpwoext*".dot"
  # fid = open(dotfile,"w")
  # write(fid,Graphs.to_dot(fgd.g))
  # close(fid)
  DFG.toDotFile(fgl, dotfile)
  show ? (@async run(`$(engine) $(dotfile) -T$(fext) -o $(filepath)`)) : nothing

  try
    viewerapp != nothing ? (@async run(`$(viewerapp) $(filepath)`)) : nothing
  catch e
    @warn "not able to show $(filepath) with viewerapp=$(viewerapp). Exception e=$(e)"
  end
  nothing
end


"""
    $(SIGNATURES)

Return the last up message stored in `cliq` of Bayes (Junction) tree.
"""
function upMsg(cliq::Graphs.ExVertex)
  @warn "deprecated upMsg, use getUpMsg instead"
  getData(cliq).upMsg
end
function upMsg(btl::BayesTree, sym::Symbol)
  @warn "deprecated upMsg, use getUpMsg instead"
  upMsg(whichCliq(btl, sym))
end

"""
    $(SIGNATURES)

Return the last down message stored in `cliq` of Bayes (Junction) tree.
"""
function dwnMsg(cliq::Graphs.ExVertex)
  @warn "deprecated dwnMsg, use getDwnMsgs instead"
  getData(cliq).dwnMsg
end
function dwnMsg(btl::BayesTree, sym::Symbol)
  @warn "deprecated dwnMsg, use getDwnMsgs instead"
  dwnMsg(whichCliq(btl, sym))
end



"""
$(SIGNATURES)

Add a node (variable) to a graph. Use this over the other dispatches.

DEPRECATED: use addVarialbe! instead.
"""
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::T;
                  N::Int=100,
                  autoinit::Bool=true,  # does init need to be separate from ready? TODO
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  uid::Int=-1,
                  smalldata=""  ) where {T <:InferenceVariable}
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!( fg,
                       lbl,
                       softtype,
                       N=N,
                       autoinit=autoinit,  # does init need to be separate from ready? TODO
                       ready=ready,
                       dontmargin=dontmargin,
                       labels=labels,
                       uid=uid,
                       smalldata=smalldata )
end
function addNode!(fg::FactorGraph,
                  lbl::Symbol,
                  softtype::Type{<:InferenceVariable};
                  N::Int=100,
                  autoinit::Bool=true,
                  ready::Int=1,
                  dontmargin::Bool=false,
                  labels::Vector{<:AbstractString}=String[],
                  uid::Int=-1,
                  # dims::Int=-1,
                  smalldata=""  )
  #
  @warn "IIF.addNode!(..) is being deprecated, use IIF.addVariable!(..) instead."
  return addVariable!(fg,
                      lbl,
                      softtype,
                      N=N,
                      autoinit=autoinit,
                      ready=ready,
                      dontmargin=dontmargin,
                      labels=labels,
                      uid=uid,
                      smalldata=smalldata  )
end

"""
    $(SIGNATURES)

Return all elements `ls(fg)` as tuples, or nodes connected to the a specific element, eg. `ls(fg, :x1)
"""
function ls(fgl::FactorGraph, lbl::Symbol; ring::Int=1)
  @warn "Deprecated, please use DFG.ls"
  # TODO ring functionality must still be implemented
  lsa = Symbol[]
  # v = nothing
  if haskey(fgl.IDs, lbl)
    id = fgl.IDs[lbl]
  else
    return lsa
  end
  # this is unnecessary
  v = getVert(fgl,id, api=api)
  for outn in api.outneighbors(fgl, v)
    # if outn.attributes["ready"] = 1 && outn.attributes["backendset"]=1
      push!(lsa, Symbol(outn.label))
    # end
  end
  return lsa
end
ls(fgl::FactorGraph, lbl::T) where {T <: AbstractString} = ls(fgl, Symbol(lbl))

"""
    $(SIGNATURES)

Experimental union of elements version of ls(::FactorGraph, ::Symbol).  Not mean't to replace broadcasting `ls.(fg, [:x1;:x2])`
"""
function ls(fgl::FactorGraph,
            lbls::Vector{Symbol};
            ring::Int=1)
  @warn "Deprecated, please use DFG.ls"
  union(ls.(fgl, lbls, ring=ring, api=api)[:]...)
end

"""
    $(SIGNATURES)

List the nodes in a factor graph.

# Examples
```julia-repl
ls(fg)
```
"""
function ls(fgl::FactorGraph; key1='x', key2='l')
  @warn "Deprecated, please use DFG.ls"
  k = collect(keys(fgl.IDs))
  x, l = String[], String[]
  xval, lval = Int[], Int[]
  xstr, lstr = String[], String[]
  xvalnested, lvalnested = String[], String[]
  xstrnested, lstrnested = String[], String[]
  canparse1, canparse2 = true,true
  nestedparse1, nestedparse2 = true, true
  idx = 0
  for id in k
    idx += 1
    idstr = string(id)
    # val = parse(Int,kstr[2:end]) # TODO: handle non-int labels
    node_idx = idstr[2:end]
    canparse = allnums(node_idx)
    nested = isnestednum(node_idx)
    if idstr[1] == key1
      keystr = string(key1,node_idx)
      if canparse
        push!(xstr, keystr)
        push!(xval, parse(Int, node_idx))
      elseif nested
        push!(xvalnested, node_idx)
        push!(xstrnested, string(node_idx))
      else
        push!(x,keystr)
      end
    elseif idstr[1] == key2
      keystr = string(key2,node_idx)
      if canparse
        push!(lstr, keystr)
        push!(lval, parse(Int, node_idx))
      elseif nested
        push!(lstrnested, keystr)
        push!(lvalnested, string(node_idx))
      else
        push!(l,string(key2,node_idx))
      end
    end
  end
  x1 = xstr[sortperm(xval)]
  x2 = xstrnested[sortnestedperm(xvalnested)]
  x = [x1; x2; sort(x)]

  l1 = lstr[sortperm(lval)]
  l2 = lstrnested[sortnestedperm(lvalnested)]
  l = [l1; l2; sort(l)]

  xx = Symbol.(x)
  ll = Symbol.(l)
  return xx, ll #return poses, landmarks
end

lsf(fgl::FactorGraph) = collect(keys(fgl.fIDs))

"""
    $(SIGNATURES)

List factors in a factor graph.

# Examples
```julia-repl
lsf(fg, :x1)
```
"""
function lsf(fgl::FactorGraph, lbl::Symbol)
  @warn "Deprecated, please use DFG.lsf"
  lsa = Symbol[]
  if haskey(fgl.fIDs, lbl)
    id = fgl.fIDs[lbl]
  else
    return lsa
  end
  v = getVert(fgl, id, api=api) # fgl.g.vertices[id] #fgl.f[id]
  for outn in api.outneighbors(fgl, v) # out_neighbors(v, fgl.g)
    push!(lsa, Symbol(outn.label))
  end
  return lsa
end


"""
    $(SIGNATURES)

List factors in a factor graph.

# Examples
```julia-repl
lsf(fg)
```
"""
lsf(fgl::FactorGraph, lbl::T) where {T <: AbstractString} = lsf(fgl,Symbol(lbl))

function lsf(fgl::FactorGraph,
             mt::Type{T}) where {T <: FunctorInferenceType}
  @warn "Deprecated, please use DFG.lsf"
  syms = Symbol[]
  for (fsym,fid) in fgl.fIDs
    if typeof(getfnctype(fgl, fid, api=api))==T
      push!(syms, fsym)
    end
  end
  return syms
end

function lsf(fgl::FactorGraph)
  @warn "Deprecated, please use DFG.lsf"
  collect(keys(fgl.fIDs))
end

"""
    $SIGNATURES

List vertices two neighbors deep.
"""
function ls2(fgl::FactorGraph, vsym::Symbol)
  @warn "Deprecated, please use DFG.ls2"
  xxf = ls(fgl, vsym)
  xlxl = Symbol[]
  for xf in xxf
    xx = lsf(fgl,xf)
    xlxl = union(xlxl, xx)
  end
  xlxl = setdiff(xlxl, [vsym])
  return xlxl
end

function initializeNode!(fgl::G,
                         sym::Symbol;
                         N::Int=100  ) where G <: AbstractDFG
  #
  @warn "initializeNode! has been deprecated in favor of initVariable!"
  initVariable!(fgl,sym,N=N )
end
