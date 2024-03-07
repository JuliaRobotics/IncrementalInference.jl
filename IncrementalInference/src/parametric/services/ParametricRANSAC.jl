
#TODO upstream to DFG
function listNeighborhood(fg, variableFactorLabels, distance; filterOrphans=true)
    allvarfacs = getNeighborhood(fg, variableFactorLabels, distance)
    variableLabels = intersect(listVariables(fg), allvarfacs)
    factorLabels = intersect(listFactors(fg), allvarfacs)
    if filterOrphans
        filter!(factorLabels) do lbl
            issubset(getVariableOrder(fg, lbl), variableLabels)
        end
    end
    return  variableLabels, factorLabels
end

function calcResidualInliers(subfg, faclabels; kσ=1)
    varlabels =  setdiff(getNeighborhood(subfg, faclabels, 1), faclabels)

    vars = getVariable.(subfg, varlabels)
    M, varTypes, vartypeslist = IIF.buildGraphSolveManifold(vars)
    varIntLabel, varlabelsAP = IIF.getVarIntLabelMap(vartypeslist)

    p0 = map(varlabelsAP) do label
        getVal(subfg, label, solveKey = :parametric)[1]
    end

    # create an ArrayPartition{CalcFactorResidual} for faclabels
    calcfacs = IIF.CalcFactorResidualAP(subfg, faclabels, varIntLabel)

    # remember res = cf.sqrt_iΣ * factor res with calcfacs, so should be sigma scaled
    res = reduce(vcat, map(f -> f(p0), Vector(calcfacs)))
    # findfirst(==(median(res)), res)

    res_labels = getproperty.(Vector(calcfacs), :faclbl)
    # findfirst(==(:cD1l2000cD1l2000_cD1l5611f1), res_labels)

    # 3
    inlierlabels = deepcopy(res_labels)
    for (i,l) in (enumerate(faclabels))
        if abs(res[i]) > kσ
            setdiff!(inlierlabels, [l])
        end
    end
    
    return res_labels.=>res, inlierlabels
end

function solveGraphParametricRansac!(
    fg,
    comb_on_fac_labels,
    min_factors = 3,
    include_vars = Symbol[];
    n_iters = 50,
    stop_ratio = 0.7, 
    kwargs...
)
    # intersect!(comb_on_fac_labels, listNeighbors(fg, ls(fg, r"^x")[1]))
    all_fac_combinations = combinations(comb_on_fac_labels, min_factors)

    do_combinations = length(all_fac_combinations) < n_iters ?
        collect(all_fac_combinations) : rand(collect(all_fac_combinations), n_iters)
        # collect(all_fac_combinations) : all_fac_combinations[rand(1:length(all_fac_combinations), n_iters)]

    best_ratio = 0.0
    best_inlierlabels = Symbol[]
    for (i, maybe_inliers) = enumerate(do_combinations)

        solveVariableLabels, solveFactorLabels = listNeighborhood(fg, union(include_vars, maybe_inliers), 2)
        
        # TODO do better
        subfg = deepcopy(fg)
        # @info solveFactorLabels
        # M, v, r, Σ = IIF.solve_RLM(fg, variableLabels, factorLabels;
        try
            IIF.solveGraphParametric!(subfg, solveVariableLabels, solveFactorLabels; kwargs...)
            # IIF.solveGraphParametric!(subfg, kwargs...)
        catch er
            @warn "Error on iter $i" er
            continue
        end

        residuals, inlierlabels = calcResidualInliers(subfg, comb_on_fac_labels)

        ratio_inliers = length(inlierlabels) ./ length(comb_on_fac_labels)

        if ratio_inliers > best_ratio
            best_ratio = ratio_inliers
            @info "new best $best_ratio"
            best_inlierlabels = inlierlabels
        end
        if ratio_inliers > stop_ratio
            @info "stop ratio met $best_ratio"
            break
        end
        # res_vals = last.(residuals)

    end
    try 
        solveVariableLabels, solveFactorLabels = listNeighborhood(fg, union(include_vars, best_inlierlabels), 2)
        IIF.solveGraphParametric!(fg, solveVariableLabels, solveFactorLabels; kwargs...)
    catch er
        @error "solveGraphParametric of inliers failed" er
    end

    return best_inlierlabels
end

if false
## get factor residuals to solve with
comb_on_fac_labels = lsfTypesDict(subfg)[:Pose2Point2Bearing]
# intersect!(faclabels, listNeighbors(subfg, ls(subfg, r"^cD1l\d+$")[1]))
stopping_criterion=StopAfterIteration(100) | StopWhenGradientNormLess(1e-12) | StopWhenStepsizeLess(1e-12)

inlierlabels = solveGraphParametricRansac!(subfg, comb_on_fac_labels;
    n_iters = 500,
    stopping_criterion,
    # debug,
    is_sparse=false,
    damping_term_min=1e-12,
    finiteDiffCovariance=true
)

end