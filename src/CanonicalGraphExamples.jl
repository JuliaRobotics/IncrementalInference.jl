
export loadCanonicalFG_Kaess, loadCanonicalFG_TestSymbolic


"""
    $SIGNATURES

Canonical example from literature, Kaess, et al.: ISAM2, IJRR, 2011.

Notes
- Paper variable ordering: p = [:l1;:l2;:x1;:x2;:x3]
"""
function loadCanonicalFG_Kaess()
  fg = initfg()

  addVariable!(fg,:x1, ContinuousScalar)
  addFactor!(fg, [:x1;], Prior(Normal()))


  addVariable!(fg,:x2, ContinuousScalar)
  addFactor!(fg,[:x1, :x2], LinearConditional(Normal()))

  addVariable!(fg, :x3, ContinuousScalar)
  addFactor!(fg,[:x2,:x3],LinearConditional(Normal()))

  addVariable!(fg, :l1, ContinuousScalar)
  addFactor!(fg, [:x1,:l1], LinearConditional(Normal()) )
  addFactor!(fg, [:x2,:l1], LinearConditional(Normal()) )

  addVariable!(fg, :l2, ContinuousScalar)
  addFactor!(fg, [:x3,:l2], LinearConditional(Normal()))

  return fg
end


"""
    $SIGNATURES

Canonical example from literature, Kaess, et al.: ISAM2, IJRR, 2011.

Notes
- Paper variable ordering: p = [:l1;:l2;:x1;:x2;:x3]
"""
function loadCanonicalFG_TestSymbolic()
  fg = initfg()

  addVariable!(fg, :x1, ContinuousScalar)
  addVariable!(fg, :x2, ContinuousScalar)
  addVariable!(fg, :x3, ContinuousScalar)
  addVariable!(fg, :x4, ContinuousScalar)
  addVariable!(fg, :x5, ContinuousScalar)
  addVariable!(fg, :l1, ContinuousScalar)
  addVariable!(fg, :l2, ContinuousScalar)
  addVariable!(fg, :l3, ContinuousScalar)

  addFactor!(fg, [:x1;:l1], LinearConditional(Normal()))
  addFactor!(fg, [:x1;:x2], LinearConditional(Normal()))
  addFactor!(fg, [:x2;:l1], LinearConditional(Normal()))
  addFactor!(fg, [:x2;:x3], LinearConditional(Normal()))
  addFactor!(fg, [:x3;:x4], LinearConditional(Normal()))
  addFactor!(fg, [:x4;:l2], LinearConditional(Normal()))
  addFactor!(fg, [:x4;:x5], LinearConditional(Normal()))
  addFactor!(fg, [:l2;:x5], LinearConditional(Normal()))
  addFactor!(fg, [:x4;:l3], LinearConditional(Normal()))
  addFactor!(fg, [:x5;:l3], LinearConditional(Normal()))

  return fg
end
