using AMD
using IncrementalInference
using Test

@testset "Test mcs, rcm, and mmd orderings" begin

fg = generateGraph_Kaess(graphinit=false)

vo = getEliminationOrder(fg, ordering=:mcs)
@test length(vo) == length(Set(vo))
@test length(vo) == length(ls(fg))

vo = getEliminationOrder(fg, ordering=:rcm)
@test length(vo) == length(Set(vo))
@test length(vo) == length(ls(fg))

vo = getEliminationOrder(fg, ordering=:mmd)
@test length(vo) == length(Set(vo))
@test length(vo) == length(ls(fg))

end
