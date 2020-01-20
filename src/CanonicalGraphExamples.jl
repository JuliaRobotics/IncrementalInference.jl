
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
  addFactor!(fg, [:x1;], Prior(Normal()), autoinit=graphinit)


  addVariable!(fg,:x2, ContinuousScalar)
  addFactor!(fg,[:x1, :x2], LinearConditional(Normal()), autoinit=graphinit)

  addVariable!(fg, :x3, ContinuousScalar)
  addFactor!(fg,[:x2,:x3],LinearConditional(Normal()), autoinit=graphinit)

  addVariable!(fg, :l1, ContinuousScalar)
  addFactor!(fg, [:x1,:l1], LinearConditional(Normal()) , autoinit=graphinit)
  addFactor!(fg, [:x2,:l1], LinearConditional(Normal()) , autoinit=graphinit)

  addVariable!(fg, :l2, ContinuousScalar)
  addFactor!(fg, [:x3,:l2], LinearConditional(Normal()), autoinit=graphinit)

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

  addFactor!(fg, [:x1;:l1], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x1;:x2], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x2;:l1], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x2;:x3], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x3;:x4], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x4;:l2], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x4;:x5], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:l2;:x5], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x4;:l3], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x5;:l3], LinearConditional(Normal()), autoinit=graphinit)

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

  addFactor!(fg, [:x0], Prior(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x0;:x1], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x1;:x2], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x2;:x3], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x3;:x4], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x4;:x5], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x5;:x6], LinearConditional(Normal()), autoinit=graphinit)

  addVariable!(fg, :l1, ContinuousScalar)
  addFactor!(fg, [:x0;:l1], LinearConditional(Normal()), autoinit=graphinit)
  addFactor!(fg, [:x6;:l1], LinearConditional(Normal()), autoinit=graphinit)

  return fg
end



"""
    $SIGNATURES
Continuous, linear scalar and multivariate test graph generation. Follows a line
with the pose id equal to the ground truth.
"""
function generateCanonicalFG_lineStep(lineLength::Int;
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
                    params=SolverParams())
                    # params=SolverParams(algorithms=[:default, :parametric]))

    vtype = (vardims == 1) ? ContinuousScalar() : ContinuousMultivariate(vardims)

    fg = LightDFG{SolverParams}( params=params)
    # fg = GraphsDFG{SolverParams}( params=params)

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
            addVariable!(fg, Symbol("x",i), vtype, autoinit = graphinit)
            (i in posePriorsAt) && addFactor!(fg, [Symbol("x",i)], Prior(xNoise(i, σ_pose_prior)))
            # "odo" type
            (i > 0) && addFactor!(fg, [Symbol("x",i-poseEvery); Symbol("x",i)], LinearConditional(xNoise(poseEvery, σ_pose_pose)))
        end


        if landmarkEvery != 0 && mod(i,landmarkEvery) == 0
            push!(lm, i)
            addVariable!(fg, Symbol("lm",i), vtype, autoinit = graphinit)
            (i in landmarkPriorsAt) && addFactor!(fg, [Symbol("lm",i)], Prior(xNoise(i, σ_lm_prior)))
        end
    end

    #add landmarks sightings
    for xi in x, lmi in lm
        dist = lmi - xi
        if abs(dist) < sightDistance
            # @info "adding landmark lm$lmi to x$xi with dist $dist"
            addFactor!(fg, [Symbol("x",xi); Symbol("lm",lmi)], LinearConditional(xNoise(dist, σ_pose_lm)))
        end
    end

    return fg
end
