

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
