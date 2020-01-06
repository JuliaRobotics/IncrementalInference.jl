# NOTE extra file to allow for esier merging later.

export lineStepFG

"""
    $SIGNATURES
Continuous, linear scalar and multivariate test graphs generation. Follows a line
with the pose id equal to the ground truth.
"""
function lineStepFG(lineLength::Int;
                    poseEvery::Int=2,
                    landmarkEvery::Int=4,
                    posePriorsAt = Int[0],
                    landmarkPriorsAt = Int[],
                    sightDistance::Int=4,
                    vardims = 1,
                    noisy = false,
                    graphinit = false,
                    σ_pose_prior = 1.0,
                    σ_lm_prior = 1.0,
                    σ_pose_pose = 1.0,
                    σ_pose_lm = 1.0,
                    params=SolverParams(algorithms=[:default, :parametric]))

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
