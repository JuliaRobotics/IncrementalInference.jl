
export loadCanonicalFG_Kaess, loadCanonicalFG_TestSymbolic, loadCanonicalFG_CaesarRing1D


"""
    $SIGNATURES

Canonical example from literature, Kaess, et al.: ISAM2, IJRR, 2011.

Notes
- Paper variable ordering: p = [:l1;:l2;:x1;:x2;:x3]
"""
function loadCanonicalFG_Kaess(;graphinit::Bool=false)
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
function loadCanonicalFG_TestSymbolic(;graphinit::Bool=false)
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
function loadCanonicalFG_CaesarRing1D(;graphinit::Bool=false)

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
