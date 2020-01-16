using IncrementalInference
using Test


@testset "Test ccolamd for constrained variable ordering" begin

fg = generateCanonicalFG_Kaess(graphinit=false)

vo = getEliminationOrder(fg, constraints=[:x3], ordering=:ccolamd)

@test vo[end] == :x3
@test length(vo) == length(ls(fg))

vo = getEliminationOrder(fg, constraints=[:l2], ordering=:ccolamd)

@test vo[end] == :l2


vo = getEliminationOrder(fg, constraints=[:x3;:l2], ordering=:ccolamd)

@test intersect(vo[end-1:end], [:x3;:l2]) |> length == 2


end
