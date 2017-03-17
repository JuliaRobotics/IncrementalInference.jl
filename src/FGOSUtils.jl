# Factor Graph OS type utilities


function convertfunctionnode(fgl::FactorGraph, fsym::Symbol, api=localapi)
  fid = fgl.fIDs[fsym]
  fnc = getfnctype(fgl, fid)
  tfn = typeof(fnc)
  nms = split(string(tfn),'.')
  typn = nms[end]
  modu = nms[1] != typn ? nms[1] : nothing  # careful with the module name
  pnm = string("Packed",typn)
  usrtyp = eval(parse(pnm))
  cfnd = convert(PackedFunctionNodeData{usrtyp},getData(fgl, fid, api=api) )
  # pfnc = convert(usrtyp, fnc)
  return cfnd, usrtyp
end

"""
    encodefg(fgl::FactorGraph)

Make a full memory copy of the graph and encode all complicated function node
types with assumed to exist convert to 'Packed<type>' formats. Same converters
as used for database persistence storage with CloudGraphs.jl.
"""
function encodefg(fgl::FactorGraph; api::DataLayerAPI=localapi)
  fgs = deepcopy(fgl)
  fgs.cg = nothing
  fgs.registeredModuleFunctions = nothing
  fgs.g = Graphs.incdict(Graphs.ExVertex,is_directed=false)
  for (vsym,vid) in fgl.IDs
    cpvert = deepcopy(getVert(fgl, vid, api=api))
    api.addvertex!(fgs, cpvert) #, labels=vnlbls)  # currently losing labels
  end

  for (fsym,fid) in fgs.fIDs
    data,ftyp = convertfunctionnode(fgl, fsym)
    # data = FunctionNodeData{ftyp}(Int64[], false, false, Int64[], m, gwpf)
    newvert = ExVertex(fid,string(fsym))
    for (key,val) in getVert(fgl,fid,api=api).attributes
      newvert.attributes[key] = val
    end
    setData!(newvert, data)
    api.addvertex!(fgs, newvert)
  end
  fgs.g.inclist = typeof(fgl.g.inclist)()

  # iterated over all edges
  for (eid, edges) in fgl.g.inclist
    fgs.g.inclist[eid] = Vector{typeof(edges[1])}()
    for ed in edges
      newed = Graphs.Edge(ed.index,
          fgs.g.vertices[ed.source.index],
          fgs.g.vertices[ed.target.index]  )
      push!(fgs.g.inclist[eid], newed)
    end
  end

  return fgs
end


"""
    savefgjld(fgl::FactorGraph; file::AbstractString="tempfg.jld")

Save mostly complete Factor Graph type by converting complicated FunctionNodeData
types to 'Packed' types using user supplied converters.
"""
function savefgjld(fgl::FactorGraph;
      file::AbstractString="tempfg.jld")
  fgs = encodefg(fgl)
  @save file fgs
  return file
end












#
