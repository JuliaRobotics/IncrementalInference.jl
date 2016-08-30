using ProtoBuf


type MyType
  x::Array{Float64,1}
  d::Dict{Int64,Int64}
  MyType() = new()
  MyType(a::Array{Float64,1},b::Dict{Int64,Int64}) = new(a,b)
end

type CompType
  mt::MyType
  str::ASCIIString
  CompType() = new()
  CompType(a::MyType,s::ASCIIString) = new(a,s)
end

iob  = PipeBuffer()

mt = MyType(rand(3),Dict{Int,Int}(1 => 10, 2 => 20))
s = ASCIIString("Hello world.")
ct = CompType(mt,s)

writeproto(iob, ct)
dd = readproto(iob, CompType())

@time writeproto(iob, ct)
@time dd = readproto(iob, CompType())
