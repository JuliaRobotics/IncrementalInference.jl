

function helloworld()
  "Hello world"
end

function atest!(rr::Dict{Int,Future})
  rr[1] = remotecall(helloworld, 1)
  nothing
end

RR = Dict{Int, Future}()
@sync begin
  @async atest!(RR)
end

@show fetch(RR[1])
