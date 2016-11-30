

function helloworld(val)
  @show "Hello world $(val)"
end

function atest!(rr::Dict{Int,Future}, val::Int)
  println("val is $(val)")
  sleep(1.0)
  rr[val] = remotecall(helloworld, 1, val)
  nothing
end

RR = Dict{Int, Future}()
iter = 1
@sync begin
  while iter <= 2
    @async atest!(RR, iter)
    iter+=1
  end
end

@show fetch(RR[1])
@show fetch(RR[2])
