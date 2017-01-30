

# DX = [tx,ty]
function odoAdd(X::Array{Float64,1}, DX::Array{Float64,1})
    # TODO -- this really needs to be cleaned up and optimized! Rather do separate Point2Point2 constraint.
    len = length(X)
    A = eye(3)
    #@show A[2,3] = X[2]
    A[1:len,3] = X

    B = eye(3)
    B[1:len,3] = DX
    #B[2,3] = DX[2]
    retval = zeros(len,1)
    AB = A*B
    retval[1:len,1] = AB[1:len,3]
    #retval[2,1] = AB[2,3]
    return retval
end

function evalPotential(odom::Odo, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int=100)
    warn("evalPotential(::Odo..) deprecated, use functor instead")
    rz,cz = size(odom.Zij)
    # implicit equation portion -- bi-directional pairwise function
    if Xid == Xi[1].index #odom.
        #Z = (odom.Zij\eye(rz)) # this will be used for group operations
        Z = - odom.Zij
        Xval = getVal(Xi[2])
    elseif Xid == Xi[2].index
        Z = odom.Zij
        Xval = getVal(Xi[1])
    else
        error("Bad evalPairwise Odo")
    end

    r,c = size(Xval)
    RES = zeros(r,c*cz)
    # increases the number of particles based on the number of modes in the measurement Z
    for i in 1:(c*cz)
        ent = diag(odom.Cov).*randn(size(vec(Z[:,floor(Int,i/(c+1)+1)])))
        # RES[:,i] = Xval[:,i] + ent+Z[:,floor(Int,i/(c+1)+1)]
        RES[:,i] = odoAdd(Xval[:,i], ent+Z[:,floor(Int,i/(c+1)+1)])
    end
    return RES
end



function evalPotential(odom::OdoMM, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int=100)
    rz,cz = size(odom.Zij)
    # implicit equation portion -- bi-directional pairwise function
    if Xid == Xi[1].index #odom.
        #Z = (odom.Zij\eye(rz)) # this will be used for group operations
        Z = -odom.Zij

        XvalMM = Array{Array{Float64,2},1}(2)
        for i in 2:3
            XvalMM[i-1] = getVal(Xi[i])
        end

        len1 = size(XvalMM[1],2)
        len2 = size(XvalMM[2],2)
        if len1 != len2
            error("odoMM -- must be the same size $(len1), $(len2)")
        #    tmp = kde!(XvalMM[2][:,:], "lcv")
        #    XvalMM[2], = sample(p,len1)
        #elseif len2 > len1
        #    tmp = kde!(XvalMM[1][:,:], "lcv")
        #    XvalMM[1], = sample(p,len2)
        end

        len = len1
        s = rand(len).>0.5 # assume fixed weight for now
        p = (1:len)[s]
        np = (1:len)[!s]
        Xval = zeros(size(XvalMM[1]))
        Xval[:,p] = XvalMM[1][:,p]
        Xval[:,np] = XvalMM[2][:,np]
    elseif Xid == Xi[2].index || Xid == Xi[3].index
        Z = odom.Zij
        Xval = getVal(Xi[1])
    else
        error("Bad evalPairwise OdoMM")
    end

    r,c = size(Xval)
    RES = zeros(r,c*cz)
    # increases the number of particles based on the number of modes in the measurement Z
    for i in 1:(c*cz)
        ent = diag(odom.Cov).*randn(size(vec(Z[:,floor(Int,i/(c+1)+1)])))
        RES[:,i] = odoAdd(Xval[:,i], ent+Z[:,floor(Int,i/(c+1)+1)])
    end
    return RES
end


function rangeAdd(X::Array{Float64,1}, DX::Array{Float64,1})
    A = eye(2)
    A[1,2] = X[1]
    B = eye(2)
    B[1,2] = DX[1]
    return [(A*B)[1,2]]
end

function evalPotential(rang::Ranged, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int=100)
    if Xid == Xi[1].index #rang.
        Z = -rang.Zij
        Xval = getVal(Xi[2])
    elseif Xid == Xi[2].index
        Z = rang.Zij
        Xval = getVal(Xi[1])
    else
        error("Bad evalPairwise Ranged")
    end

    r,c = size(Xval)
    cz = size(rang.Zij,1)
    RES = zeros(r,c*cz)

    for i in 1:(c*cz) # for each mode in the measurement
        ent = rang.Cov[1]*randn(size(vec(Z[:,floor(Int,i/(c+1)+1)])))
        RES[:,i] = rangeAdd(Xval[:,i], ent+Z[floor(Int,i/(c+1)+1)])
    end
    return RES
end


function getSample(obs::Obsv2, N::Int=1)
  pd = kde!(obs.pts, obs.bws[:,1])
  return (KernelDensityEstimate.sample(pd,N)[1],)
end
# TODO -- this may be obsolete, investigate further and remove
function evalPotential(obs::Obsv2, Xi::Array{Graphs.ExVertex,1}; N::Int64=100)#, from::Int64)
  return getSample(obs, N)
end

function evalPotentialSpecific(fnc::Function, Xi::Array{Graphs.ExVertex,1}, typ::Singleton, solvefor::Int64; N::Int64=100)
  outpts = fnc(typ, Xi, N=N) # , solvefor
  # outpts = evalPotential(typ, Xi, N=N) # , solvefor
  return outpts
end

function evalPotentialSpecific(fnc::Function, Xi::Array{Graphs.ExVertex,1}, typ::Pairwise, solvefor::Int64; N::Int64=100)
  return fnc(typ, Xi, solvefor, N=N)
  # return evalPotential(typ, Xi, solvefor)
end

# function prepareparamsarray!(ARR::Array{Array{Float64,2},1},Xi::Vector{Graphs.ExVertex}, N::Int, solvefor::Int64)
#   LEN = Int64[]
#   maxlen = N
#   count = 0
#   sfidx = 0
#   for xi in Xi
#     push!(ARR, getVal(xi))
#     len = size(ARR[end], 2)
#     push!(LEN, len)
#     if len > maxlen
#       maxlen = len
#     end
#     count += 1
#     if xi.index == solvefor
#       sfidx = xi.index
#     end
#   end
#   SAMP=LEN.<maxlen
#   for i in 1:count
#     if SAMP[i]
#       ARR[i] = KernelDensityEstimate.sample(getKDE(Xi[i]), maxlen)[1]
#     end
#   end
#   # we are generating a proposal distribution, not direct replacement for existing
#   ARR[sfidx] = deepcopy(ARR[sfidx])
#   return maxlen, sfidx
# end


# currently monster of a spaghetti code mess, multiple types interacting at the same
# time. Safest at this point in development is to just get this code running and
# will refactor once the dust has settled.
# current code is the result of several unit tests between IIF and RoME.jl
# a unified test to follow -- after which refactoring can start
function evalPotentialSpecific{T <: FunctorPairwise}(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int64;
      N::Int64=100  )
  #
  # TODO -- enable partial constraints

  # TODO -- this part can be collapsed into common generic solver component, could be constructed and maintained at addFactor! time
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end
  # gwp.zDim[sfidx] ??
  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)
  # and return complete fr/gwp

  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  for gwp.particleidx in 1:maxlen
    # gwp(x, res)
    numericRootGenericRandomizedFnc!( fr )
        # r = nlsolve( gwp, ARR[sfidx][:,gwp.particleidx] )
        # remember this is a deepcopy of original sfidx, since we are generating a proposal distribution
        # and not directly replacing the existing variable belief estimate
        # gwp.params[gwp.varidx][:,gwp.particleidx] = r.zero[:]
  end

  return gwp.params[gwp.varidx]
  # return evalPotential(typ, Xi, solvefor)
end

function evalPotentialSpecific{T <: FunctorPairwiseMinimize}(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      gwp::GenericWrapParam{T},
      solvefor::Int64;
      N::Int64=100  )
  #
  # TODO -- this part can be collapsed into common generic solver component, could be constructed and maintained at addFactor! time
  ARR = Array{Array{Float64,2},1}()
  maxlen, sfidx = prepareparamsarray!(ARR, Xi, N, solvefor)
  gwp.params = ARR
  gwp.varidx = sfidx
  gwp.measurement = gwp.samplerfnc(gwp.usrfnc!, maxlen)
  zDim = size(gwp.measurement[1],1) # TODO -- zDim aspect desperately needs to be redone
  if gwp.specialzDim
    zDim = gwp.usrfnc!.zDim[sfidx]
  end

  fr = FastRootGenericWrapParam{T}(gwp.params[sfidx], zDim, gwp)

  # TODO -- once Threads.@threads have been optmized JuliaLang/julia#19967, also see area4 branch
  for gwp.particleidx in 1:maxlen
    # gwp(x, res)
    # implement minimization here
    res = zeros(fr.xDim)
    gg = (x) -> fr.gwp(x, res)
    r = optimize(  gg, fr.X[1:fr.xDim,fr.gwp.particleidx] )
    # TODO -- clearly lots of optmization to be done here
    fr.Y[1:fr.xDim] = r.minimizer
    fr.X[:,fr.gwp.particleidx] = fr.Y

    # error("not implemented yet")
  end
  return gwp.params[gwp.varidx]
end

function evalPotentialSpecific{T <: FunctorSingleton}(
      fnc::T,
      Xi::Vector{Graphs.ExVertex},
      generalwrapper::GenericWrapParam{T},
      solvefor::Int64;
      N::Int64=100  )
  #
  generalwrapper.measurement = generalwrapper.samplerfnc(generalwrapper.usrfnc!, N)
  return generalwrapper.measurement[1]
end

# Multiple dispatch occurs internally, resulting in factor graph potential evaluations
function evalFactor2(fgl::FactorGraph, fct::Graphs.ExVertex, solvefor::Int64; N::Int64=100)
  # return evalPotential(fct.attributes["data"].fnc, solvefor) #evalPotential(fct.attributes["fnc"], solvefor)

  # TODO -- this build up of Xi is excessive and should be reduced
  # could happen at addFactor time
  Xi = Graphs.ExVertex[]
  for id in fct.attributes["data"].fncargvID
    push!(Xi, dlapi.getvertex(fgl,id)) # TODO localapi
  end
  # lookup now used for getSample
  # modulefnc = fgl.registeredModuleFunctions[fct.attributes["data"].frommodule]
  fnctype = fct.attributes["data"].fnc
  return evalPotentialSpecific(fnctype.usrfnc!, Xi, fnctype, solvefor, N=N)
  # return evalPotentialSpecific(modulefnc, Xi, fnctype, solvefor, N=N)
end

function findRelatedFromPotential(fg::FactorGraph, idfct::Graphs.ExVertex, vertid::Int64, N::Int64) # vert
  # if vert.index == vertid
  # assuming it is properly initialized TODO
    ptsbw = evalFactor2(fg, idfct, vertid, N=N);
    # NOTE -- disable this validation test
    # sum(abs(ptsbw)) < 1e-14 ? error("findRelatedFromPotential -- an input is zero") : nothing

    # TODO -- better to upsample before the projection
    Ndim = size(ptsbw,1)
    Npoints = size(ptsbw,2)
    # Assume we only have large particle population sizes, thanks to addNode!
    p = kde!(ptsbw, "lcv")
    if Npoints != N # this is where we control the overall particle set size
        p = resample(p,N)
    end
    return p
  # end
  # return Union{}
end
