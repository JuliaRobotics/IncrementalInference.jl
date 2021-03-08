
export generateCanonicalFG_Kaess, generateCanonicalFG_TestSymbolic, generateCanonicalFG_CaesarRing1D, generateCanonicalFG_lineStep

"""
    $SIGNATURES

Canonical example from literature, Kaess, et al.: ISAM2, IJRR, 2011.

Notes
- Paper variable ordering: p = [:l1;:l2;:x1;:x2;:x3]
"""
function generateCanonicalFG_Kaess(;graphinit::Bool=false)
  fg = initfg()

  addVariable!(fg,:x1, ContinuousScalar)
  addFactor!(fg, [:x1;], Prior(Normal()), graphinit=graphinit)


  addVariable!(fg,:x2, ContinuousScalar)
  addFactor!(fg,[:x1, :x2], LinearRelative(Normal()), graphinit=graphinit)

  addVariable!(fg, :x3, ContinuousScalar)
  addFactor!(fg,[:x2,:x3],LinearRelative(Normal()), graphinit=graphinit)

  addVariable!(fg, :l1, ContinuousScalar)
  addFactor!(fg, [:x1,:l1], LinearRelative(Normal()) , graphinit=graphinit)
  addFactor!(fg, [:x2,:l1], LinearRelative(Normal()) , graphinit=graphinit)

  addVariable!(fg, :l2, ContinuousScalar)
  addFactor!(fg, [:x3,:l2], LinearRelative(Normal()), graphinit=graphinit)

  return fg
end


"""
    $SIGNATURES

Canonical example introduced by Borglab.

Notes
- Known variable ordering: p = [:x1; :l3; :l1; :x5; :x2; :l2; :x4; :x3]
"""
function generateCanonicalFG_TestSymbolic(;graphinit::Bool=false)
  fg = initfg()

  addVariable!(fg, :x1, ContinuousScalar)
  addVariable!(fg, :x2, ContinuousScalar)
  addVariable!(fg, :x3, ContinuousScalar)
  addVariable!(fg, :x4, ContinuousScalar)
  addVariable!(fg, :x5, ContinuousScalar)
  addVariable!(fg, :l1, ContinuousScalar)
  addVariable!(fg, :l2, ContinuousScalar)
  addVariable!(fg, :l3, ContinuousScalar)

  addFactor!(fg, [:x1;:l1], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x1;:x2], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x2;:l1], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x2;:x3], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x3;:x4], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x4;:l2], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x4;:x5], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:l2;:x5], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x4;:l3], LinearRelative(Normal()), graphinit=graphinit)
  addFactor!(fg, [:x5;:l3], LinearRelative(Normal()), graphinit=graphinit)

  return fg
end



"""
    $SIGNATURES

Canonical example introduced originally as Caesar Hex Example.

Notes
- Paper variable ordering: p = [:x0;:x2;:x4;:x6;:x1;:l1;:x5;:x3;]
"""
function generateCanonicalFG_CaesarRing1D(;graphinit::Bool=false)

    fg = initfg()

    addVariable!(fg, :x0, ContinuousScalar)
    addVariable!(fg, :x1, ContinuousScalar)
    addVariable!(fg, :x2, ContinuousScalar)
    addVariable!(fg, :x3, ContinuousScalar)
    addVariable!(fg, :x4, ContinuousScalar)
    addVariable!(fg, :x5, ContinuousScalar)
    addVariable!(fg, :x6, ContinuousScalar)

    addFactor!(fg, [:x0], Prior(Normal()), graphinit=graphinit)
    addFactor!(fg, [:x0;:x1], LinearRelative(Normal()), graphinit=graphinit)
    addFactor!(fg, [:x1;:x2], LinearRelative(Normal()), graphinit=graphinit)
    addFactor!(fg, [:x2;:x3], LinearRelative(Normal()), graphinit=graphinit)
    addFactor!(fg, [:x3;:x4], LinearRelative(Normal()), graphinit=graphinit)
    addFactor!(fg, [:x4;:x5], LinearRelative(Normal()), graphinit=graphinit)
    addFactor!(fg, [:x5;:x6], LinearRelative(Normal()), graphinit=graphinit)

    addVariable!(fg, :l1, ContinuousScalar)
    addFactor!(fg, [:x0;:l1], LinearRelative(Normal()), graphinit=graphinit)
    addFactor!(fg, [:x6;:l1], LinearRelative(Normal()), graphinit=graphinit)

    return fg
end



"""
    $SIGNATURES
Continuous, linear scalar and multivariate test graph generation. Follows a line
with the pose id equal to the ground truth.
"""
function generateCanonicalFG_lineStep(  lineLength::Int;
                                        poseEvery::Int=2,
                                        landmarkEvery::Int=4,
                                        posePriorsAt = Int[0],
                                        landmarkPriorsAt = Int[],
                                        sightDistance::Int=4,
                                        vardims = 1,
                                        noisy = false,
                                        graphinit = false,
                                        σ_pose_prior = 0.1,
                                        σ_lm_prior = 0.1,
                                        σ_pose_pose = 0.1,
                                        σ_pose_lm = 0.1,
                                        solverParams=SolverParams())
                                        # solverParams=SolverParams(algorithms=[:default, :parametric]))

    vtype = (vardims == 1) ? ContinuousScalar() : ContinuousEuclid(vardims)

    fg = LightDFG{SolverParams}( solverParams=solverParams)

    function xNoise(i::Int, σ::Float64=1.0)
        if (vardims == 1)
            return noisy ? Normal(σ*randn() + i, σ) : Normal(0.0*randn() + i, σ)
        else
            return noisy ? MvNormal(σ*randn(vardims) .+ i, σ) : MvNormal(0.0*randn(vardims) .+ i, σ)
        end
    end

    x = Int[]
    lm = Int[]
    for i=0:lineLength
        if mod(i,poseEvery) == 0
            push!(x, i)
            addVariable!(fg, Symbol("x",i), vtype) #, autoinit = graphinit)
            (i in posePriorsAt) && addFactor!(fg, [Symbol("x",i)], Prior(xNoise(i, σ_pose_prior)), graphinit=graphinit)
            # "odo" type
            (i > 0) && addFactor!(fg, [Symbol("x",i-poseEvery); Symbol("x",i)], LinearRelative(xNoise(poseEvery, σ_pose_pose)), graphinit=graphinit)
        end


        if landmarkEvery != 0 && mod(i,landmarkEvery) == 0
            push!(lm, i)
            addVariable!(fg, Symbol("lm",i), vtype) #, autoinit = graphinit)
            (i in landmarkPriorsAt) && addFactor!(fg, [Symbol("lm",i)], Prior(xNoise(i, σ_lm_prior)), graphinit=graphinit)
        end
    end

    #add landmarks sightings
    for xi in x, lmi in lm
        dist = lmi - xi
        if abs(dist) < sightDistance
            # @info "adding landmark lm$lmi to x$xi with dist $dist"
            addFactor!(fg, [Symbol("x",xi); Symbol("lm",lmi)], LinearRelative(xNoise(dist, σ_pose_lm)), graphinit=graphinit)
        end
    end

    return fg
end

export generateCanonicalFG_EuclidDistance
"""
    $SIGNATURES
Generate a EuclidDistance test graph where 1 landmark position is unknown. 
"""
function generateCanonicalFG_EuclidDistance(points::Vector{Vector{Float64}} = [[100.,0],[0.,100]];
                                            dist = 100.0, 
                                            σ_prior=1.0, σ_dist=1.0,
                                            N=100,
                                            graphinit=false)
    #
    dims = length(points[1])
    fg = initfg()
    fg.solverParams.N = N
    fg.solverParams.graphinit = graphinit


    for (i,p) in enumerate(points)
    xlbl = Symbol("x",i)
    addVariable!(fg, xlbl, ContinuousEuclid{dims})
    addFactor!(fg, [xlbl], Prior(MvNormal(p, σ_prior*ones(dims))))
    end

    addVariable!(fg, :l1, ContinuousEuclid{dims})

    for i=1:length(points)
    xlbl = Symbol("x",i)
    addFactor!(fg, [xlbl;:l1], EuclidDistance(Normal(dist, σ_dist)))
    end

    return fg
end