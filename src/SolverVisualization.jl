import KernelDensityEstimate.kde!

function drawOneMC!(cliqMC::CliqGibbsMC, minmax, mcmc=0; offs=2.0)

    if mcmc>minmax[4]
      minmax[4]=mcmc
    end

    i = 0.0
    for prod in cliqMC.prods
        prodVal = kde!(prod.product,"lcv") #cliqMC.prods[1]
        plotKDEProd!([prodVal;prod.potentials],minmax, h=-i*offs, mcmc=mcmc)
        i += 1.0
    end

end

function drawMCMCDebug(cliq; offs=2.0)
    println("$(cliq.attributes["label"])")
    cliqDbg = cliq.attributes["debug"]
    MCd = 100.0/length(cliqDbg.mcmc)
    MCo=0.0
    minmax=[99999,-99999,-10,-99999.0]
    for mc in cliqDbg.mcmc
        drawOneMC!(mc, minmax, MCo, offs=offs)
        MCo+=MCd
    end

    n = 3
    y = linspace(minmax[1], minmax[2], n)
    x = linspace(minmax[3],minmax[4],n)
    xgrid = repmat(x',n,1)
    ygrid = repmat(y,1,n)
    z = zeros(n,n)

    for i in 1:length(cliqDbg.mcmc[1].prods)
      surf(xgrid,ygrid,z-(i-1)*offs,alpha=0.04, linewidth=0.0)
    end
end


function drawTreeUpwardMsgs(fgl::FactorGraph, bt::BayesTree; N=300)
    len = length(bt.cliques)-1
    vv = Array{Gadfly.Compose.Context,1}(len)
    #r = Array{RemoteRef,1}(len)
    i = 0
    for cliq in bt.cliques
        if cliq[1] == 1 println("No upward msgs from root."); continue; end
        @show cliq[2].attributes["label"]
        i+=1
        vv[i] = drawHorDens(fgl, cliq[2].attributes["debug"].outmsg.p, N)
    end
    #[r[j] = @spawn drawCliqueMsgs(bt.cliques[j+1]) for j in 1:len]
    #[vv[j] = fetch(r[j]) for j in 1:len]
    vv
end

function drawFrontalDens(fg::FactorGraph, bt::BayesTree;
                        N=300, gt=Union{})
    len = length(bt.cliques)
    vv = Array{Gadfly.Compose.Context,1}(len)
    i = 0
    for cliq in bt.cliques
        #@show cliq[2].attributes["label"]
        lenfr = length(cliq[2].attributes["frontalIDs"])

        p = Array{BallTreeDensity,1}(lenfr)
        j=0
        #pvals = Array{Array{Float64,2},1}(lenfr)
        gtvals = Dict{Int,Array{Float64,2}}()
        lbls = ASCIIString[]

        for frid in cliq[2].attributes["frontalIDs"]
            j+=1
            p[j] = kde!(fg.v[frid].attributes["val"])
            #pvals[j] = fg.v[frid].attributes["val"]

            if gt!=Union{}
              gtvals[j] = gt[fg.v[frid].attributes["label"]]
              #push!(gtvals, gt[fg.v[frid].attributes["label"]][1])
              #push!(gtvals, gt[fg.v[frid].attributes["label"]][2])
            end
            push!(lbls, fg.v[frid].attributes["label"])

        end

        #r = Array{RemoteRef,1}(lenfr)
        #[r[j] = @spawn kde!(pvals[j]) for j in 1:lenfr]
        #[p[j] = fetch(r[j]) for j in 1:lenfr]

        i+=1
        if length(gtvals) > 0
          #gtvals = reshape(gtvals,2,round(Int,length(gtvals)/2))'
          vv[i] = drawHorDens(p, N, gt=gtvals, lbls=lbls)
        else
          vv[i] = drawHorDens(p, N,lbls=lbls)
        end
    end
    #
    vv
end

# for some reason we still need msgPlots of right size in the global for these functions to work.
# precall drawTreeUpwardMsgs or drawFrontalDens to make this work properly TODO
function vstackedDensities(msgPlots)
    #msgPlots = f(fg, bt) # drawTreeUpwardMsgs
    evalstr = ""
    for i in 1:length(msgPlots)
        evalstr = string(evalstr, ",msgPlots[$(i)]")
    end
    eval(parse(string("vstack(",evalstr[2:end],")")))
end
#
# ## TODO -- you were here with port starboard lines
# function stbPrtLineLayers!(pl, Xpp, Ypp, Thpp)
#
#     l = 5.0
#     lnstpr = [0.0;l;0.0]
#     lnstpg = [0.0;-l;0.0]
#
#     Rd  =SE2(lnstpr)
#     Gr = SE2(lnstpg)
#
#     for i in 1:length(Xpp)
#       lnstt = [Xpp[i];Ypp[i];Thpp[i]]
#       Ps = SE2(lnstt)
#       lnr = se2vee(Ps*Rd)
#       lng = se2vee(Ps*Gr)
#       xsr = [Xpp[i];lnr[1]]
#       ysr = [Ypp[i];lnr[2]]
#       xsg = [Xpp[i];lng[1]]
#       ysg = [Ypp[i];lng[2]]
#
#       push!(pl.layers, layer(x=xsr, y=ysr, Geom.path(), Gadfly.Theme(default_color=colorant"red", line_width=1.5pt))[1] )
#       push!(pl.layers, layer(x=xsg, y=ysg, Geom.path(), Gadfly.Theme(default_color=colorant"green", line_width=1.5pt))[1] )
#     end
#     nothing
# end
#
# # function lblsFromTo(from,to)
# #   lbls=ASCIIString[]
# #   [push!(lbls, "$(i)") for i in from:to]
# #   return lbls
# # end
#
# function drawPoses(fg::FactorGraph; from::Int64=0,to::Int64=99999999,
#                     meanmax=:max, lbls=true, drawhist=true)
#     #Gadfly.set_default_plot_size(20cm, 30cm)
#     Xp,Yp = get2DPoseSamples(fg, from=from, to=to)
#     Xpp = Float64[]; Ypp=Float64[]; Thpp=Float64[]; LBLS=ASCIIString[];
#     if meanmax == :mean
#       Xpp,Ypp, Thpp, LBLS = get2DPoseMeans(fg, from=from, to=to)
#     elseif meanmax == :max
#       Xpp,Ypp, Thpp, LBLS = get2DPoseMax(fg, from=from, to=to)
#     end
#
#     # lbls = lblsFromTo(1,length(Xpp))
#     psplt = Union{}
#     if lbls
#       psplt = Gadfly.plot(
#       Gadfly.layer(x=Xpp,y=Ypp,label=LBLS,Geom.path(), Theme(line_width=2pt), Geom.label)
#       )
#     else
#       psplt = Gadfly.plot(
#       Gadfly.layer(x=Xpp,y=Ypp,Geom.path(), Theme(line_width=2pt))
#       )
#     end
#     stbPrtLineLayers!(psplt, Xpp, Ypp, Thpp)
#     if drawhist
#       push!(psplt.layers,  Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)[1] )#(xbincount=100, ybincount=100))
#     end
#     psplt
# end
#
# function drawLandms(fg::FactorGraph;
#               from::Int64=0, to::Int64=99999999, minnei::Int64=0,
#               meanmax=:max,lbls=true,showmm=false,drawhist=true,c="red",MM=Union{})
#     #Gadfly.set_default_plot_size(20cm, 30cm)
#     Xp,Yp = get2DLandmSamples(fg, from=from, to=to)
#     Xpp = Float64[]; Ypp=Float64[]; Thpp=Float64[]; lblstags=ASCIIString[];
#     if meanmax==:mean
#       Xpp,Ypp, t, lbltags = get2DLandmMeans(fg, from=from, to=to)
#     elseif meanmax==:max
#       Xpp,Ypp, t, lbltags = get2DLandmMax(fg, from=from, to=to,showmm=showmm,MM=MM)
#     end
#
#     if lbls
#       psplt = Gadfly.plot(
#       Gadfly.layer(x=Xpp,y=Ypp, label=lbltags, Geom.point, Theme(line_width=2pt, default_color=parse(Colorant,c)), Geom.label)
#       # ,Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)#(xbincount=100, ybincount=100)
#       )
#     else
#       psplt = Gadfly.plot(
#       Gadfly.layer(x=Xpp,y=Ypp, Geom.point, Theme(line_width=2pt, default_color=parse(Colorant,c)))
#       )
#     end
#
#     if drawhist
#       push!(psplt.layers, Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)[1])#(xbincount=100, ybincount=100)
#     end
#
#     psplt
# end
#
# function drawPosesLandms(fg::FactorGraph;
#                     from::Int64=0, to::Int64=99999999, minnei::Int64=0,
#                     meanmax=:max,lbls=true,drawhist=true, MM=Union{})
#   p = drawPoses(fg, from=from,to=to,meanmax=meanmax,lbls=lbls,drawhist=drawhist)
#   pl = drawLandms(fg, from=from, to=to, minnei=minnei,lbls=lbls,drawhist=drawhist, MM=MM)
#   for l in pl.layers
#     push!(p.layers, l)
#   end
#   return p
# end
#
# function drawSubmaps(fgl::FactorGraph, fromto::Array{Int,2};
#                     m1hist=false,m2hist=false,m3hist=false, showmm=false)
#   p = drawLandms(fgl, from=fromto[1,1], to=fromto[1,2], drawhist=m1hist, showmm=showmm)
#   if size(fromto,1) >1
#     p2 = drawLandms(fgl, from=fromto[2,1], to=fromto[2,2], drawhist=m2hist,c="blue", showmm=showmm)
#     for l in p2.layers
#       push!(p.layers, l)
#     end
#   end
#   if size(fromto,1) >2
#     p3 = drawLandms(fgl, from=fromto[3,1], to=fromto[3,2], drawhist=m3hist,c="magenta", showmm=showmm)
#     for l in p3.layers
#       push!(p.layers, l)
#     end
#   end
#   return p
# end
#
# function drawSubmaps(fgl::FactorGraph, fromto::Array{Int,1}; spread::Int=25,
#                     m1hist=false,m2hist=false,m3hist=false, showmm=false)
#   ft = zeros(Int,length(fromto),2)
#   for i in 1:length(fromto)
#     ft[i,1] = fromto[i]-spread; ft[i,2] = fromto[i]+spread;
#   end
#   drawSubmaps(fgl, ft, m1hist=m1hist, m2hist=m2hist, m3hist=m3hist, showmm=showmm)
# end
#
# # function getKDEMax(p::BallTreeDensity;N=200)
# #   m = zeros(p.bt.dims)
# #   for i in 1:p.bt.dims
# #     mm = marginal(p,[i])
# #     rangeV = getKDERange(mm)
# #     X = linspace(rangeV[1],rangeV[2],N)
# #     yV = evaluateDualTree(mm,X)
# #     m[i] = X[findfirst(yV,maximum(yV))]
# #   end
# #   return m
# # end
#
function investigateMultidimKDE(p::BallTreeDensity, p0::BallTreeDensity)
    co = ["black"; "blue"]
    h = Union{}
    x = plotKDE([marginal(p,[1]); marginal(p0,[1])], c=co )
    y = plotKDE([marginal(p,[2]); marginal(p0,[2])], c=co )
    if p.bt.dims >= 3
      th = plotKDE([marginal(p,[3]); marginal(p0,[3])], c=co )
      h = hstack(x,y,th)
    else
      h = hstack(x,y)
    end

    return h
end


function investigateMultidimKDE(p::Array{BallTreeDensity,1})
    co = ["black"; "blue"; "green"; "red"; "magenta"; "cyan"; "cyan1"; "cyan2";
    "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta";
    "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"]
    # compute all the marginals
    Pm = Array{Array{BallTreeDensity,1},1}()
    push!(Pm,stackMarginals(p,1)) #[marginal(p[1],[1]); marginal(p[2],[1])]
    push!(Pm,stackMarginals(p,2)) #[marginal(p[1],[2]); marginal(p[2],[2])]

    h = Union{}
    x = plotKDE(Pm[1], c=co )
    y = plotKDE(Pm[2], c=co )
    if p[1].bt.dims >= 3
      #Pm3 = [marginal(p[1],[3]); marginal(p[2],[3])]
      push!(Pm,stackMarginals(p,3)) # [marginal(p[1],[3]); marginal(p[2],[3])]
      th = plotKDE(Pm[3], c=co )
      h = hstack(x,y,th)
    else
      h = hstack(x,y)
    end
    return h
end

function investigateMultidimKDE(p::BallTreeDensity)
    x = plotKDE(marginal(p,[1]) )
    y = plotKDE(marginal(p,[2]) )
    if p.bt.dims >= 3
      th = plotKDE(marginal(p,[3]) )
      return hstack(x,y,th)
    end
    return hstack(x,y)
end

function vArrPotentials(potens::Dict{Int,EasyMessage})
  vv = Array{Gadfly.Compose.Context,1}(length(potens))
  i = 0
  oned=false
  for p in potens
      i+=1
      pb = kde!(p[2].pts, p[2].bws)
      if size(p[2].pts,1) > 3
        # vv[i] = investigatePoseKDE(pb)
        error("can't handle higher dimensional plots here yet")
      elseif size(p[2].pts,1) > 1
        vv[i] = investigateMultidimKDE(pb)
      else
        vv[i] = plotKDE(pb)
      end
  end
  # if oned
    return vv
  # else
  #   println("wrong one")
  #   return vv2
  # end
end

function kde!(em::EasyMessage)
  return kde!(em.pts,em.bws)
end

function draw(em::EasyMessage;xlbl="X")
  p = Union{}
  if size(em.pts,1) == 1
    p=plotKDE(kde!(em),xlbl=xlbl)
  else
    p=investigatePoseKDE(kde!(em))
  end
  return p
end

function whosWith(cliq::Graphs.ExVertex)
  println("$(cliq.attributes["label"])")
  for pot in cliq.attributes["potentials"]
      println("$(pot[2])")
  end
  nothing
end

function whosWith(bt::BayesTree, frt::ASCIIString)
    whosWith(whichCliq(bt,frt))
end


function drawUpMsgAtCliq(fg::FactorGraph, cliq::Graphs.ExVertex)
    for id in keys(cliq.attributes["debug"].outmsg.p)
        print("$(fg.v[id].attributes["label"]), ")
    end
    println("")
    sleep(0.1)
    potens = cliq.attributes["debug"].outmsg.p
    vArrPotentials(potens)
end

function drawUpMsgAtCliq(fg::FactorGraph, bt::BayesTree, lbl::ASCIIString)
    drawUpMsgAtCliq(fg, whichCliq(bt, lbl) )
end

# function drawDwnMsgAtCliq(fg::FactorGraph, cliq::Graphs.ExVertex)
function dwnMsgsAtCliq(fg::FactorGraph, cliq::Graphs.ExVertex)
  for id in keys(cliq.attributes["debugDwn"].outmsg.p)
      print("$(fg.v[id].label), ")
  end
  println("")
  sleep(0.1)
  potens = cliq.attributes["debugDwn"].outmsg.p
  potens
end

function dwnMsgsAtCliq(fg::FactorGraph, bt::BayesTree, lbl::ASCIIString)
    dwnMsgsAtCliq(fg, whichCliq(bt, lbl) )
end

function drawPose2DMC!(plots::Array{Gadfly.Compose.Context,1}, cliqMC::CliqGibbsMC)

    for prod in cliqMC.prods
        prodVal = kde!(prod.product,"lcv") #cliqMC.prods[1]
        push!(plots, investigatePoseKDE([prodVal;prod.potentials]) )
    end
    vstackedPlots(plots)
end

function mcmcPose2D!(plots::Array{Gadfly.Compose.Context,1}, cliqDbg::DebugCliqMCMC, iter::Int=1)
    # for mc in cliqDbg.mcmc
    mc = cliqDbg.mcmc[iter]
    v = drawPose2DMC!(plots, mc)
    # end
    return v
end

function drawUpMCMCPose2D!(plots::Array{Gadfly.Compose.Context,1}, cliq::Graphs.ExVertex, iter::Int=1)
    whosWith(cliq)
    cliqDbg = cliq.attributes["debug"]
    sleep(0.1)
    mcmcPose2D!(plots, cliqDbg, iter)
end

function drawUpMCMCPose2D!(plots::Array{Gadfly.Compose.Context,1}, bt::BayesTree, frt::ASCIIString, iter::Int=1)
    drawUpMCMCPose2D!(plots, whichCliq(bt,frt), iter)
end

function drawDwnMCMCPose2D!(plots::Array{Gadfly.Compose.Context,1}, cliq::Graphs.ExVertex, iter::Int=1)
    whosWith(cliq)
    cliqDbg = cliq.attributes["debugDwn"]
    sleep(0.1)
    mcmcPose2D!(plots, cliqDbg, iter)
end

function drawDwnMCMCPose2D!(plots::Array{Gadfly.Compose.Context,1}, bt::BayesTree, frt::ASCIIString, iter::Int=1)
    drawDwnMCMCPose2D!(plots, whichCliq(bt,frt), iter)
end

function drawLbl(fg::FactorGraph, lbl::ASCIIString)
    v = fg.v[fg.IDs[lbl]]
    investigatePoseKDE(kde!(getVal(v)))
end

function predCurrFactorBeliefs(fgl::FactorGraph, fc::Graphs.ExVertex)
  # TODO update to use ls and lsv functions
  prjcurvals = Dict{ASCIIString, Array{BallTreeDensity,1}}()
  for v in out_neighbors(fc, fgl.g)
    pred = kde!(evalFactor2(fgl, fc, v.index))
    curr = kde!(getVal(v))
    prjcurvals[v.attributes["label"]] = [curr; pred]
  end
  return prjcurvals, collect(keys(prjcurvals))
end


function drawHorDens(fgl::FactorGraph, pDens::Dict{Int,EasyMessage}, N=200)
  p = BallTreeDensity[]
  lbls = ASCIIString[]
  for pd in pDens
    push!(p, kde!(pd[2].pts,pd[2].bws))
    push!(lbls, fgl.v[pd[1]].attributes["label"])
  end
  @show lbls
  drawHorDens(p,N,lbls=lbls)
end

function drawHorBeliefsList(fgl::FactorGraph, lbls::Array{ASCIIString,1};
                        nhor::Int=-1,gt=Union{},N::Int=200, extend=0.1)
  len = length(lbls)
  pDens = BallTreeDensity[]
  for lb in lbls
    pt = getVal(fgl.v[fgl.IDs[lb]])
    push!(pDens, kde!(pt,"lcv"))
  end

  if nhor<1
    nhor = round(Int,sqrt(len))
  end
  vlen = ceil(Int, len/nhor)
  vv = Array{Gadfly.Compose.Context,1}(vlen)
  conslb = deepcopy(lbls)
  vidx = 0
  for i in 1:nhor:len
    pH = BallTreeDensity[]
    gtvals = Dict{Int,Array{Float64,2}}()
    labels = ASCIIString[]
    for j in 0:(nhor-1)
      if i+j <= len
        push!(pH, pDens[i+j])
        push!(labels, lbls[i+j])
        if gt != Union{} gtvals[j+1] = gt[lbls[i+j]]  end
      end
    end
    vidx+=1
    if gt !=Union{}
      vv[vidx] = drawHorDens(pH, N=N, gt=gtvals, lbls=labels, extend=extend)
    else
      vv[vidx] = drawHorDens(pH, N=N, lbls=labels, extend=extend)
    end
  end
  vv
end

function drawFactorBeliefs(fgl::FactorGraph, flbl::ASCIIString)
  if !haskey(fgl.fIDs, flbl)
    println("no key $(flbl)")
    return nothing
  end
  # for fc in fgl.f
    # if fc[2].attributes["label"] == flbl
    fc = fgl.f[fgl.fIDs[flbl]]
      prjcurvals, lbls = predCurrFactorBeliefs(fgl, fc)
      if length(lbls) == 3
        return vstack(
        investigatePoseKDE(prjcurvals[lbls[1]]),
        investigatePoseKDE(prjcurvals[lbls[2]]),
        investigatePoseKDE(prjcurvals[lbls[3]]),
        )
      elseif length(lbls) == 2
        return vstack(
        investigatePoseKDE(prjcurvals[lbls[1]]),
        investigatePoseKDE(prjcurvals[lbls[2]]),
        )
      elseif length(lbls) == 1
        return investigatePoseKDE(prjcurvals[lbls[1]])
      end

    # end
  # end
  nothing
end

function localProduct(fgl::FactorGraph, lbl::ASCIIString;N=300)
  arr = Array{BallTreeDensity,1}()
  cf = ls(fgl, lbl)
  pp = Union{}
  for f in cf
    push!(arr, kde!(evalFactor2(fgl, fgl.f[fgl.fIDs[f]], fgl.IDs[lbl])))
  end
  if length(arr)>1
    Ndims = arr[1].bt.dims
    dummy = kde!(rand(Ndims,N),[1.0]);
    pGM, = prodAppxMSGibbsS(dummy, arr, Union{}, Union{}, 40)
    pp = kde!(pGM,"lcv")
  end
  return pp,arr
end

function drawLocalProduct(fgl::FactorGraph, lbl::ASCIIString;N=300)
  arr = Array{BallTreeDensity,1}()
  push!(arr, getVertKDE(fgl, lbl)) #kde!(getVal(fgl.v[fgl.IDs[lbl]]))
  pp, parr = localProduct(fgl, lbl, N=N)
  if pp != Union{}
    push!(arr,pp)
    for a in parr
      push!(arr, a)
    end
  end
  return investigatePoseKDE(arr)
end

function saveplot(pl;name="pl",frt=:png,w=25cm,h=25cm,nw=false,fill=true)
  if frt==:png
    Gadfly.draw(PNG(string(name,".png"),w,h),pl)
    if fill run(`composite $(name).png plB.png $(name).png`) end
    if !nw run(`eog $(name).png`) end
  end
  if frt==:pdf
    Gadfly.draw(PDF(string(name,".pdf"),w,h),pl)
    if !nw run(`evince $(name).pdf`) end
  end
  nothing
end

function animateVertexBelief(FGL::Array{FactorGraph,1}, lbl;nw=false)
  len = length(FGL)
  [saveplot(drawLocalProduct(FG[i],lbl),h=15cm,w=30cm,name="gifs/pl$(i)",nw=true) for i=1:len];
  run(`convert -delay 100 gifs/pl*.png result.gif`)
  if !nw run(`eog result.gif`) end
  nothing
end

function ls(fgl::FactorGraph, lbl::ASCIIString)
  ls = ASCIIString[]
  v = Union{}
  if haskey(fgl.IDs, lbl)
    id = fgl.IDs[lbl]
  else
    return ls
  end
  v = fgl.v[id]
  for outn in out_neighbors(v, fgl.g)
    push!(ls, outn.label)
  end
  return ls
end

function ls(fgl::FactorGraph)
  k = collect(keys(fgl.IDs))
  l = Int[]
  x = Int[]
  for kk in k
    val = parse(Int,kk[2:end])
    if kk[1] == 'l'
      push!(l,val)
    elseif kk[1] == 'x'
      push!(x,val)
    end
  end
  l = sort(l)
  x = sort(x)
  ll = Array{ASCIIString,1}(length(l))
  xx = Array{ASCIIString,1}(length(x))
  for i in 1:length(l)
    ll[i] = string("l",l[i])
  end
  for i in 1:length(x)
    xx[i] = string("x",x[i])
  end
  return xx,ll
end

function lsv(fgl::FactorGraph, lbl::ASCIIString)
  ls = ASCIIString[]
  v = Union{}
  if haskey(fgl.fIDs, lbl)
    id = fgl.fIDs[lbl]
  else
    return ls
  end
  v = fgl.f[id]
  for outn in out_neighbors(v, fgl.g)
    push!(ls, outn.label)
  end
  return ls
end



function fixRotWrapErr!(RT::Array{Float64,1})

  for i in 1:length(RT)
    if RT[i] > pi
      RT[i] = abs(RT[i]-2.0*pi)
    end
  end
  nothing
end

function asyncUniComp(fgl::FactorGraph, isamdict::Dict{Int,Array{Float64,1}})
  r,rt,lb = computeGraphResiduals(fgl,isamdict);
  fixRotWrapErr!(rt)
  return [sqrt(Base.mean(r.^2));maximum(abs(r));sqrt(Base.mean(rt.^2));maximum(rt)]
end

function unimodalCompare(FGL::Array{FactorGraph,1},isamdict::Dict{Int,Array{Float64,1}})
  len = length(FGL)
  RMS = Float64[]
  MAX = Float64[]
  RMSth = Float64[]
  MAXth = Float64[]

  rr = RemoteRef[]

  for fgl in FGL
    push!(rr, remotecall(uppA(),asyncUniComp, fgl, isamdict))
  end

  for r in rr
    err = fetch(r)
    push!(RMS, err[1])
    push!(MAX, err[2])
    push!(RMSth, err[3])
    push!(MAXth, err[4])
  end

  x=0:(len-1)
  df1 = DataFrame(x=x, y=RMS, label="rms")
  df2 = DataFrame(x=x, y=MAX, label="max")
  df3 = DataFrame(x=x, y=RMSth*180.0/pi, label="rmsth")
  df4 = DataFrame(x=x, y=MAXth*180.0/pi, label="maxth")
  df = vcat(df1, df2)
  dfth = vcat(df3,df4)

  return df,dfth
end

function asyncAnalyzeSolution(fgl::FactorGraph, lbl::ASCIIString)
  pp, arr = localProduct(fgl, lbl)
  lpm = getKDEMax(pp)
  em = getKDEMax(getVertKDE(fgl,lbl))
  err1 = norm(lpm[1:2]-em[1:2])
  err2 = 0.0
  if lbl[1]=='x'
    err2 = abs(lpm[3]-em[3])
  end
  return [err1;err2]
end

function analyzeSolution(FGL::Array{FactorGraph,1},fggt=Union{})
  len = length(FGL)
  RMS = Float64[]
  MAX = Float64[]
  RMSth = Float64[]
  MAXth = Float64[]
  for fgl in FGL
    xLB, lLB = ls(fgl)
    ERR = Float64[]
    ERRth = Float64[]
    ALB = [xLB;lLB]
    rr = RemoteRef[]
    for lbl in ALB
      push!(rr, remotecall(uppA(),asyncAnalyzeSolution, fgl, lbl))
      # err
      # push!(ERR, err[1])
      # if lbl[1]=='x'
      #   push!(ERRth, err[2])
      # end
    end

    idx = 1
    for r in rr
      err = fetch(r)
      push!(ERR, err[1])
      if ALB[idx][1]=='x'
        push!(ERRth, err[2])
      end
      idx += 1
    end
    push!(RMS, sqrt(Base.mean(ERR.^2)))
    push!(MAX, maximum(abs(ERR)))
    push!(RMSth, sqrt(Base.mean(ERRth.^2)))
    push!(MAXth, maximum(ERRth))
  end

  x=0:(len-1)
  df1 = DataFrame(x=x, y=RMS, label="rms")
  df2 = DataFrame(x=x, y=MAX, label="max")
  df3 = DataFrame(x=x, y=RMSth*180.0/pi, label="rmsth")
  df4 = DataFrame(x=x, y=MAXth*180.0/pi, label="maxth")
  df = vcat(df1, df2)
  dfth = vcat(df3,df4)
  # return Gadfly.plot(
  # layer(,y=RMS,Geom.path(),Geom.point, Gadfly.Theme(default_color=parse(Colorant,"red"))),
  # layer(x=1:len,y=MAX,Geom.path(),Geom.point, Gadfly.Theme(default_color=parse(Colorant,"black"))),
  # layer(x=1:len,y=RMSth*180.0/pi,Geom.path(),Geom.point, Gadfly.Theme(default_color=parse(Colorant,"blue"))),
  # layer(x=1:len,y=MAXth*180.0/pi,Geom.path(),Geom.point, Gadfly.Theme(default_color=parse(Colorant,"green")),
  # )
  return df,dfth
end
# discrete_color_manual(colors...; levels=nothing,order=nothing) is deprecated, use color_discrete_manual(colors...; levels=levels,order=order) instead.

function drawAnalysis(df,dfth)
  return vstack(
  Gadfly.plot(df, x="x", y="y", color="label", Geom.line,
         Scale.discrete_color_manual("red","black")),
  Gadfly.plot(dfth, x="x", y="y", color="label", Geom.line,
        Scale.discrete_color_manual("red","black"))
        )
end

function getAllFGsKDEs(fgD::Array{FactorGraph,1}, vertid::Int64)
  ret = Array{BallTreeDensity,1}()
  for i in 1:length(fgD)
    V = fgD[i].v[vertid]
    push!(ret, kde!(getVal(V)))
  end
  return ret
end

function drawAllPose2DBeliefs(plots::Array{Gadfly.Compose.Context,1}, fgD::Array{FactorGraph,1})
    ids = sort(collect(keys(fgD[1].v)))
    co = ["black"; "blue"; "green"; "red"; "magenta"; "cyan"; "cyan1"; "cyan2"]
    println(co[1:length(fgD)])
    # V = fg.v[7]
    for i in ids
        # V = fgD[1].v[i]
        @show fgD[1].v[i].attributes["label"]
        # if length(fgD) >= 2
            # V0 = fgD[2].v[i]
            kdes = getAllFGsKDEs(fgD, i)
            push!(plots, investigatePoseKDE(  kdes  )) # [kde!(getVal(V)); kde!(getVal(V0))]
        # else
        #     push!(plots, investigatePoseKDE(  [kde!(getVal(V))]  ))
        # end
    end
    vstackedPlots(plots)
end

# legacy function -- use the array version instead
function drawAllPose2DBeliefs(plots::Array{Gadfly.Compose.Context,1}, fgD::FactorGraph, fgD0=Union{})
  println("WARNING: drawAllPose2DBeliefs -- legacy function -- use the array version instead.")
  if fgD0 != Union{}
    drawAllPose2DBeliefs(plots, [fgD;fgD0])
  else
    drawAllPose2DBeliefs(plots, [fgD])
  end
end

function drawComicStripLM(fgD::Array{FactorGraph,1})
    comicA = Array{Gadfly.Plot,1}()
    for fgd in fgD
        cv = drawPosesLandms(fgd)
        # cv = drawPoses(fgd)
        push!(comicA,cv)
    end
    hstack(comicA)
end

# function drawComicStripLM(fgD::Array{FactorGraph,1})
#   len = length(fgD)
#   comicA = Array{Gadfly.Plot,1}(len)
#   r = Array{RemoteRef,1}(len)
#   # for fgd in fgD
#   #     cv = drawPosesLandms(fgd)
#   #     # cv = drawPoses(fgd)
#   #     push!(comicA,cv)
#   # end
#   [r[i] = @spawn drawPosesLandms(fgD[i]) for i in 1:len]
#   [comicA[i] = fetch(r[i]) for i in 1:len]
#   hstack(comicA)
# end

function drawComicStrip(fgD::Array{FactorGraph,1})
    comicA = Array{Gadfly.Plot,1}()
    for fgd in fgD
        cv = drawPoses(fgd)
        push!(comicA,cv)
    end
    hstack(comicA)
end


function compositeComic(fnc::Function, fgGT, fgA::Array{FactorGraph,1})
    v = Union{}
    @show length(fgA)
    if length(fgA) == 2
        Gadfly.set_default_plot_size(25cm, 10cm)
        v = fnc([fgA[1:2];fgGT])
    elseif length(fgA) == 3
        Gadfly.set_default_plot_size(25cm, 20cm)
        v = vstack(fnc(fgA[1:2])
        ,fnc([fgA[3];fgGT])    )
    elseif length(fgA) == 4
        Gadfly.set_default_plot_size(25cm, 20cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc([fgA[4];fgGT])    )
    elseif length(fgA) == 7
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc([fgA[7];fgGT])    )
    elseif length(fgA) == 10
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])
        ,fnc([fgA[10];fgGT])    )
    elseif length(fgA) == 13
        Gadfly.set_default_plot_size(25cm, 30cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])
        ,fnc(fgA[10:12])
        ,fnc([fgA[13];fgGT])    )
    end
    v
end


function compositeComic(fnc::Function, fgA::Array{FactorGraph,1})
    v = Union{}
    @show length(fgA)
    if length(fgA) == 2
        Gadfly.set_default_plot_size(25cm, 10cm)
        v = fnc(fgA[1:2])
    elseif length(fgA) == 3
        Gadfly.set_default_plot_size(25cm, 20cm)
        v = fnc(fgA[1:3])
    elseif length(fgA) == 4
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:2])
        ,fnc(fgA[3:4]))
    elseif length(fgA) == 6
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])  )
    elseif length(fgA) == 9
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])    )
    elseif length(fgA) == 12
        Gadfly.set_default_plot_size(25cm, 30cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])
        ,fnc(fgA[10:12])    )
    end
    v
end
