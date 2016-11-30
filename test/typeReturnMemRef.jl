


type MyType
   arr::Array{Float64,1}
   var::Int
   str::String
   MyType(x...) = new(x[1],x[2],x[3])
   MyType() = new()
end

m = MyType(rand(3),2,"hello")

function f(M::MyType)
   mm = MyType()
   mm.arr = M.arr
   return mm
end

mmm = f(m)

mmm.arr[1] = 1.0

@test m.arr[1] == 1.0
