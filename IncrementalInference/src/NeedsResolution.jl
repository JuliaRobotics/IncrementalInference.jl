
# FIXME move to DFG
getPointDefault(V::StateType) = getPointIdentity(V)

function compare(c1::Channel, c2::Channel; skip::Vector{Symbol} = [])
  #
  TP = true
  TP = TP && c1.state == c2.state
  TP = TP && c1.sz_max == c2.sz_max
  TP = TP && c1.data |> length == c2.data |> length
  # exit early if tests already failed
  !TP && (return false)
  # now check contents of data
  for i = 1:length(c1.data)
    TP = TP && c1.data[i] == c2.data[i]
  end
  return TP
end

compare(a::Int, b::Int) = a == b
compare(a::Bool, b::Bool) = a == b
compare(a::Dict, b::Dict) = a == b
