
@info "Defining spyCliqMat(..) for visualizing association matrix of a clique in the Bayes (Junction) tree"

exportimg(pl) = Gadfly.PNG(pl)

export spyCliqMat

"""
    $SIGNATURES

Draw the clique association matrix, with keyword arguments for more or less console print outs.

Notes
* Columns are variables, rows are factors.
* Drawn from up message passing perspective.
* Blue color implies no factor association.
* Frontal, separator, and upmessages are all drawn at different intensity of red.
* Downward messages not shown, as they would just be singletons of the full separator set.
"""
function spyCliqMat(cliq::TreeClique; showmsg=true, suppressprint::Bool=false)
  mat = deepcopy(getCliqMat(cliq, showmsg=showmsg))
  # TODO -- add improved visualization here, iter vs skip
  mat = map(Float64, mat)*2.0.-1.0
  numlcl = size(getCliqAssocMat(cliq),1)
  mat[(numlcl+1):end,:] *= 0.9
  mat[(numlcl+1):end,:] .-= 0.1
  numfrtl1 = floor(Int,length(getCliqueData(cliq).frontalIDs) + 1)
  mat[:,numfrtl1:end] *= 0.9
  mat[:,numfrtl1:end] .-= 0.1
  if !suppressprint
    @show getCliqueData(cliq).itervarIDs
    @show getCliqueData(cliq).directvarIDs
    @show getCliqueData(cliq).msgskipIDs
    @show getCliqueData(cliq).directFrtlMsgIDs
    @show getCliqueData(cliq).directPriorMsgIDs
  end
  if size(mat,1) == 1
    mat = [mat; -ones(size(mat,2))']
  end
  sp = Gadfly.spy(mat)
  push!(sp.guides, Gadfly.Guide.title("$(getLabel(cliq)) || $(getCliqueData(cliq).frontalIDs) :$(getCliqueData(cliq).separatorIDs)"))
  push!(sp.guides, Gadfly.Guide.xlabel("fmcmcs $(getCliqueData(cliq).itervarIDs)"))
  push!(sp.guides, Gadfly.Guide.ylabel("lcl=$(numlcl) || msg=$(size(getCliqMsgMat(cliq),1))" ))
  return sp
end
function spyCliqMat(bt::AbstractBayesTree, lbl::Symbol; showmsg=true, suppressprint::Bool=false)
  spyCliqMat(getClique(bt,lbl), showmsg=showmsg, suppressprint=suppressprint)
end