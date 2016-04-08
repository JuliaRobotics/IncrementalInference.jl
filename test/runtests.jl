using IncrementalInference
using Base.Test

global pass=false

try
  fg = emptyFactorGraph();
  # add stuff to fg
  # tree = prepBatchTree!(fg, ordering=:qr);
  global pass=true
catch e
  global pass=false
  rethrow(e)
end

@test pass
