# test deconv functions

using IncrementalInference


fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addFactor!(fg, [:x0], Prior(Normal()))

addVariable!(fg, :x1, ContinuousScalar)
addFactor!(fg, [:x0; :x1], LinearConditional(Normal(10, 1)))


ensureAllInitialized!(fg)






#
