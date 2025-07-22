##

# using Revise
using Test
using IncrementalInference

##


@testset "Test differential factors for MKD sampling types (useMsgLikelihoods" begin
##

fg = generateGraph_CaesarRing1D()
getSolverParams(fg).useMsgLikelihoods = true

## test getSample

fct = fg[:x0x1f1]
fT = getFactorType(fct)
@test fT isa LinearRelative
@test fT.Z isa Normal

##

M = getManifold(fT)
X = sampleTangent(M, fT.Z)
@test X isa AbstractVector{<:Real}

z = sampleFactor(fct)[1]
@test z isa AbstractVector{<:Real}

##

initAll!(fg)

##

# fl = lsf(fg) |> sortDFG
X_ = approxDeconvBelief(fg,:x0f1)
X_ = approxDeconvBelief(fg,:x0x1f1)


##

eliminationOrder = [:x3,:x5,:l1,:x1,:x6,:x4,:x2,:x0]

tree = buildTreeReset!(fg, eliminationOrder)

##

# cfg = buildCliqSubgraph(fg, tree[6])
# st = IIF._generateMsgJointRelativesPriors(cfg, :default, tree[6])

cfg = buildCliqSubgraph(fg, tree[5])
st = IIF._generateMsgJointRelativesPriors(cfg, :default, tree[5])
beliefMsg5 = IIF.prepCliqueMsgUp(cfg, tree[5], :default, IIF.UPSOLVED)

cfg = buildCliqSubgraph(fg, tree[4])
st = IIF._generateMsgJointRelativesPriors(cfg, :default, tree[4])
beliefMsg4 = IIF.prepCliqueMsgUp(cfg, tree[4], :default, IIF.UPSOLVED)

# cfg = buildCliqSubgraph(fg, tree[3])
# st = IIF._generateMsgJointRelativesPriors(cfg, :default, tree[3])

cfg = buildCliqSubgraph(fg, tree[2])
@test 3 === length(ls(cfg))
@test 0 === length(lsf(cfg))

IIF.addMsgFactors!(cfg, beliefMsg4, IIF.UpwardPass)
IIF.addMsgFactors!(cfg, beliefMsg5, IIF.UpwardPass)
@test 2 === length(lsf(cfg))

##

fct = cfg[:x0x6f1]
fT = getFactorType(fct)
@test fT isa LinearRelative
@test fT.Z isa MKD

##

M = getManifold(fT.Z)
X = sampleTangent(M, fT.Z)
@test X isa AbstractVector{<:Real}

z = sampleFactor(fct)[1]
@test z isa AbstractVector{<:Real}

##


childmsgs = LikelihoodMessage[]
retdict = IIF.upGibbsCliqueDensity(cfg, tree[2], :default, childmsgs)

# st = IIF._generateMsgJointRelativesPriors(cfg, :default, tree[2])

# cfg = buildCliqSubgraph(fg, tree[1])
# st = IIF._generateMsgJointRelativesPriors(cfg, :default, tree[1])


##


getSolverParams(fg).downsolve = false
solveGraph!(fg)



##
end