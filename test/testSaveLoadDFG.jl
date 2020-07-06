
using IncrementalInference
using Test

@testset "saving to and loading from FileDFG" begin

  fg = generateCanonicalFG_Kaess()
  addVariable!(fg, :x4, ContinuousScalar)
  addFactor!(fg, [:x2;:x3;:x4], LinearConditional(Normal()), multihypo=[1.0;0.6;0.4])

  saveFolder = "/tmp/dfg_test"
  saveDFG(fg, saveFolder)
  # VERSION above 1.0.x hack required since Julia 1.0 does not seem to havfunction `splitpath`
  if v"1.1" <= VERSION
    retDFG = initfg()
    retDFG = loadDFG(saveFolder, IncrementalInference, retDFG)
    Base.rm(saveFolder*".tar.gz")

    @test symdiff(ls(fg), ls(retDFG)) == []
    @test symdiff(lsf(fg), lsf(retDFG)) == []

    @show getFactor(fg, :x2x3x4f1).solverData.multihypo
    @show getFactor(retDFG, :x2x3x4f1).solverData.multihypo

    # check for match
    @test getFactor(fg, :x2x3x4f1).solverData.multihypo - getFactor(retDFG, :x2x3x4f1).solverData.multihypo |> norm < 1e-10
    @test getFactor(fg, :x2x3x4f1).solverData.certainhypo - getFactor(retDFG, :x2x3x4f1).solverData.certainhypo |> norm < 1e-10
  end

end


#
