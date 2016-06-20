#
#
# # define the pose group
# type Odo <: Pairwise
#     #Xi::Array{Graphs.ExVertex,1}
#     Zij::Array{Float64,2} # 0rotations, 1translation in each column
#     Cov::Array{Float64,2}
#     W::Array{Float64,1}
# end
# type PackedOdo
#     #Xi::Array{Graphs.ExVertex,1}
#     vecZij::Array{Float64,1} # 0rotations, 1translation in each column
#     dimz::Int64
#     vecCov::Array{Float64,1}
#     dimc::Int64
#     W::Array{Float64,1}
# end
# function convert(::Type{Odo}, d::PackedOdo)
#   Zij = reshapeVec2Mat(d.vecZij,d.dimz)
#   Cov = reshapeVec2Mat(d.vecCov, d.dimc)
#   return Odo(Zij, Cov, d.W)
# end
# function convert(::Type{PackedOdo}, d::Odo)
#   v1 = d.Zij[:];
#   v2 = d.Cov[:];
#   return PackedOdo(v1,length(v1),
#                     v2,length(v2),
#                     d.W)
# end
# function convert(::Type{FunctionNodeData{PackedOdo}}, d::FunctionNodeData{Odo})
#   return FunctionNodeData{PackedOdo}(d.fncargvID, d.eliminated, d.potentialused,
#           convert(PackedOdo, d.fnc))
# end
# function convert(::Type{FunctionNodeData{Odo}}, d::FunctionNodeData{PackedOdo})
#   return FunctionNodeData{Odo}(d.fncargvID, d.eliminated, d.potentialused,
#           convert(Odo, d.fnc))
# end
#
#
# type OdoMM <: Pairwise
#     #Xi::Array{Graphs.ExVertex,1} # modes are Xi 2:end
#     Zij::Array{Float64,2} # 0rotations, 1translation in each column
#     Cov::Array{Float64,2}
#     W::Array{Float64,1}
# end
#
#
# type Ranged <: Pairwise
#     #Xi::Array{Graphs.ExVertex,1}
#     Zij::Array{Float64,1}
#     Cov::Array{Float64,1}
#     W::Array{Float64,1}
# end
#
# type GenericMarginal <: Pairwise
#     #Xi::Array{Graphs.ExVertex,1}
#     Zij::Array{Float64,1}
#     Cov::Array{Float64,1}
#     W::Array{Float64,1}
#     GenericMarginal() = new()
#     GenericMarginal(a,b,c,d) = new(a,b,c,d)
# end


# DX = [tx,ty]
function odoAdd(X::Array{Float64,1}, DX::Array{Float64,1})
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

function evalPotential(odom::Odo, Xi::Array{Graphs.ExVertex,1}, Xid::Int64)
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
        RES[:,i] = odoAdd(Xval[:,i], ent+Z[:,floor(Int,i/(c+1)+1)])
    end
    return RES
end

function evalPotential(odom::OdoMM, Xi::Array{Graphs.ExVertex,1}, Xid::Int64)
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

function evalPotential(rang::Ranged, Xi::Array{Graphs.ExVertex,1}, Xid::Int64)
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

#
# type Obsv <: Singleton
#     #Xi::Array{Graphs.ExVertex,1}
#     Zi::Array{Float64,1}
#     Cov::Array{Float64,1}
#     W::Array{Float64,1}
# end

# this is not right, should return many points
# function evalPotential(obs::Obsv, Xi::Array{Graphs.ExVertex,1}; N::Int64=100) #, from::Int64)
#     return obs.Cov[1]*randn()+obs.Zi
# end
#
# type Obsv2 <: Singleton
#     #Xi::Array{Graphs.ExVertex,1}
#     pts::Array{Float64,2}
#     bws::Array{Float64,2}
#     W::Array{Float64,1}
# end
# type PackedObsv2
#     #Xi::Array{Graphs.ExVertex,1}
#     vecZij::Array{Float64,1} # 0rotations, 1translation in each column
#     dimz::Int64
#     vecCov::Array{Float64,1}
#     dimc::Int64
#     W::Array{Float64,1}
# end
# function convert(::Type{Obsv2}, d::PackedObsv2)
#   Zij = reshapeVec2Mat(d.vecZij,d.dimz)
#   Cov = reshapeVec2Mat(d.vecCov, d.dimc)
#   return Obsv2(Zij, Cov, d.W)
# end
# function convert(::Type{PackedObsv2}, d::Obsv2)
#   v1 = d.pts[:];
#   v2 = d.bws[:];
#   return PackedObsv2(v1,length(v1),
#                     v2,length(v2),
#                     d.W)
# end
# function convert(::Type{FunctionNodeData{PackedObsv2}}, d::FunctionNodeData{Obsv2})
#   return FunctionNodeData{PackedObsv2}(d.fncargvID, d.eliminated, d.potentialused,
#           convert(PackedObsv2, d.fnc))
# end
# function convert(::Type{FunctionNodeData{Obsv2}}, d::FunctionNodeData{PackedObsv2})
#   return FunctionNodeData{Obsv2}(d.fncargvID, d.eliminated, d.potentialused,
#           convert(Obsv2, d.fnc))
# end


function evalPotential(obs::Obsv2, Xi::Array{Graphs.ExVertex,1}; N::Int64=300)#, from::Int64)
    # @show obs.bws, typeof(obs.bws)
    pd = kde!(obs.pts, obs.bws[:,1])
    return KernelDensityEstimate.sample(pd,N)[1]
    # return obs.Cov[1]*randn()+obs.Zi
end


function evalPotentialSpecific(Xi::Array{Graphs.ExVertex,1}, typ::Singleton, solvefor::Int64)
  return evalPotential(typ, Xi) # , solvefor
end

function evalPotentialSpecific(Xi::Array{Graphs.ExVertex,1}, typ::Pairwise, solvefor::Int64)
  return evalPotential(typ, Xi, solvefor)
end

function evalFactor2(fgl::FactorGraph, fct::Graphs.ExVertex, solvefor::Int64)
  # return evalPotential(fct.attributes["data"].fnc, solvefor) #evalPotential(fct.attributes["fnc"], solvefor)
  Xi = Graphs.ExVertex[]
  for id in fct.attributes["data"].fncargvID
    push!(Xi,fgl.v[id])
  end
  return evalPotentialSpecific(Xi, fct.attributes["data"].fnc, solvefor) #evalPotential(fct.attributes["fnc"], solvefor)
end

function findRelatedFromPotential(fg::FactorGraph, idfct::Graphs.ExVertex, vertid::Int64, N::Int64) # vert
    # if vert.index == vertid
        ptsbw = evalFactor2(fg, idfct, vertid); # idfct[2] # assuming it is properly initialized TODO
        sum(abs(ptsbw)) < 1e-14 ? error("findRelatedFromPotential -- an input is zero") : nothing

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
