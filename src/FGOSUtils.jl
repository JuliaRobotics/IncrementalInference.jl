# Factor Graph OS type utilities
#  IIF methods should direclty detect extended types from user import
# of convert in their namespace



manikde!(pts::AbstractArray{Float64,2}, vartype::InferenceVariable) = manikde!(pts, getManifolds(vartype))
manikde!(pts::AbstractArray{Float64,2}, vartype::Type{<:InferenceVariable}) = manikde!(pts, getManifolds(vartype))


"""
    $SIGNATURES

Return N=100 measurement samples for a factor in `<:AbstractDFG`.
"""
function getMeasurements(dfg::AbstractDFG, fsym::Symbol, N::Int=100)
  fnc = getFactorFunction(dfg, fsym)
  getSample(fnc, N)
end

"""
    $SIGNATURES

Get graph node (variable or factor) dimension.
"""
getDimension(var::DFGVariable) = getSofttype(var).dims
getDimension(fct::DFGFactor) = solverData(fct).fnc.zDim


"""
    $(SIGNATURES)

Return array of all variable nodes connected to the last `n` many poses (`:x*`).

Example:

```julia
# Shallow copy the tail end of poses from a factor graph `fg1`
vars = lsRear(fg1, 5)
fg1_r5 = subgraphFromVerts(fg1, vars)
```
"""
function lsRear(fgl::FactorGraph, n::Int=1)
  @warn "lsRear in current form is not gauranteed to work right, use sort(ls(fg, r\"x\")[end]) instead.  Also see sortVarNested"
  lasts = ls(fgl)[1][(end-n):end]
  syms = ls(fgl, lasts)
  union(lsf.(fgl, syms)[:]...)
end


"""
    $SIGNATURES

Return `Vector{Symbol}` of landmarks attached to vertex vsym in `fgl::FactorGraph`.
"""
function landmarks(fgl::FactorGraph, vsym::Symbol)
  fsyms = ls(fgl, vsym)
  lms = Symbol[]
  for fs in fsyms
    for varv = lsf(fgl, fs)
      if string(varv)[1] == 'l'
        push!(lms, varv)
      end
    end
  end
  lms
end



function evalLikelihood(fg::FactorGraph, sym::Symbol, point::Vector{Float64})
  p = getVertKDE(fg, sym)
  Ndim(p) == length(point) ? nothing : error("point (dim=$(length(point))) must have same dimension as belief (dim=$(Ndim(p)))")
  evaluateDualTree(p, reshape(point,:,1))[1]
end

# Evaluate the likelihood of an Array{2} of points on the marginal belief of some variable
# note the dimensions must match
function evalLikelihood(fg::FactorGraph, sym::Symbol, points::Array{Float64,2})
  p = getVertKDE(fg, sym)
  Ndim(p) == size(points,1) ? nothing : error("points (dim=$(size(points,1))) must have same dimension as belief (dim=$(Ndim(p)))")
  evaluateDualTree(p, (points))
end



function setThreadModel!(fgl::FactorGraph;model=IncrementalInference.SingleThreaded)
  for (key, id) in fgl.fIDs
    solverData(getFactor(fgl, key)).fnc.threadmodel = model
  end
  nothing
end

# not sure if and where this is still being used
function _evalType(pt::String)::Type
    try
        getfield(Main, Symbol(pt))
    catch ex
        io = IOBuffer()
        showerror(io, ex, catch_backtrace())
        err = String(take!(io))
        error("_evalType: Unable to locate factor/distribution type '$pt' in main context (e.g. do a using on all your factor libraries). Please check that this factor type is loaded into main. Stack trace = $err")
    end
end

"""
    $(SIGNATURES)

Print the maximum point values form all variables approximate marginals in the factor graph.
The full marginal can be recovered for example `X0 = getVertKDE(fg, :x0)`.
"""
function printgraphmax(fgl::FactorGraph)
    verts = union(ls(fgl)...)
    map(v -> println("$v : $(getKDEMax(getVertKDE(fgl, v)))"), verts);
end


"""
    $SIGNATURES

Return interger index of desired variable element.

Example
-------
```julia
pp = RoME.Point2()
getIdx(pp, :posY) # = 2
```

Internal Notes
--------------
- uses number i < 100 for index number, and
- uses +100 offsets to track the minibatch number of the requested dimension
"""
function getIdx(pp::T, sym::Symbol, i::Int=0)::Tuple{Int, Int} where {T <: Tuple}
  # i > 99 ? error("stop") : nothing
  i-=100
  for p in pp
    i,j = getIdx(p, sym, i)
    if i > 0
      return i, j
    end
  end
  return i,-1
end
getIdx(pp::Symbol, sym::Symbol, i::Int=0)::Tuple{Int, Int} = pp==sym ? (abs(i)%100+1, div(abs(i)-100,100)) : (i-1, div(abs(i)-100,100))
function getIdx(pp::V, sym::Symbol, i::Int=0)::Tuple{Int, Int} where {V <: InferenceVariable}
  return getIdx(pp.dimtype, sym)
end


"""
   $SIGNATURES

Display the content of `VariableNodeData` to console for a given factor graph and variable tag`::Symbol`.

Dev Notes
- TODO split as two show macros between AMP and DFG
"""
function showVariable(fgl::G, vsym::Symbol) where G <: AbstractDFG
  vert = DFG.getVariable(fg, vsym)
  vnd = solverData(vert)
  println("label: $(vert.label), exVertexId: $(vert.index)")
  println("tags: $( haskey(vert.attributes, string(:tags)) ? vert.attributes[string(:tags)] : string(:none))")
  println("size marginal samples $(size(getVal(vnd)))")
  println("kde bandwidths: $(getBW(vnd)[:,1])")
  println("kde mean: $(round.(getKDEMean(getKDE(vnd)),digits=4))")
  println("kde max: $(round.(getKDEMax(getKDE(vnd)),digits=4))")
  println()
  vnd
end

"""
    $SIGNATURES

Return `::Bool` on whether this variable has been marginalized.
"""
isMarginalized(vert::DFGVariable) = solverData(vert).ismargin
isMarginalized(dfg::AbstractDFG, sym::Symbol) = isMarginalized(DFG.getVariable(dfg, sym))



"""
    $SIGNATURES

Free all variables from marginalization.
"""
function dontMarginalizeVariablesAll!(fgl::G) where G <: AbstractDFG
  fgl.solverParams.isfixedlag = false
  fgl.solverParams.qfl = 9999999999
  for sym in ls(fgl)
    solverData(getVariable(fgl, sym)).ismargin = false
  end
  nothing
end

"""
    $SIGNATURES

Free all variables from marginalization.

Related

dontMarginalizeVariablesAll!
"""
function unfreezeVariablesAll!(fgl::G) where G <: AbstractDFG
  dontMarginalizeVariablesAll!(fgl)
end

"""
    $SIGNATURES

Reset initialization flag on all variables in `::FactorGraphs`.

Notes
- Numerical values remain, but inference will overwrite since init flags are now `false`.
"""
function resetVariableAllInitializations!(fgl::FactorGraph)
  vsyms = ls(fgl)
  for sym in vsyms
    setVariableInitialized!(getVariable(fgl, sym), :false)
  end
  nothing
end






function convert(::Type{Tuple{BallTreeDensity,Float64}},
                 p::EasyMessage )
  (AMP.manikde!(p.pts, p.bws, p.manifolds), p.inferdim)
end

function convert(::Type{EasyMessage},
                 bel::Tuple{BallTreeDensity,Float64},
                 manifolds::T) where {T <: Tuple}
  EasyMessage(getPoints(bel[1]), getBW(bel[1])[:,1], manifolds, bel[2])
end


#
