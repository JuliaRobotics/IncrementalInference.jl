
println("[TEST] Protobuf quick check on converting types and packing data")

module MyMod

  export MyType, MyTemplType, thing

  type MyType
    x::Array{Float64,1}
    val::Int64
    str::ASCIIString
    MyType() = new()
    MyType(x...) = new(x[1],x[2],x[3])
  end

  type MyTemplType{T}
    x::Array{Float64,1}
    val::Int64
    str::ASCIIString
    a::T
    MyTemplType() = new()
    MyTemplType(x...) = new(x[1],x[2],x[3],x[4])
  end

  function thing()
    println("Thing is running")
    nothing
  end

end

module ProtoTest
  using ProtoBuf
  export protostring, nothingspecialprotoread, registerType, runOutsideFunction, specialprotoread

  regTypes = Dict{ASCIIString, Type}()

  function protostring(datas)
    iob = PipeBuffer()
    @show typeof(datas)
    writeproto(iob, datas)
    return iob
  end

  function nothingspecialprotoread(iob, instance)
    return readproto(iob, instance)
  end

  typeConterters = Dict{AbstractString,Function}()
  function registerType(str::ASCIIString, spType::Type; converter::Union{Function,Union}=Union{})
    ProtoTest.regTypes[str] = spType
    if converter!=Union{}
      typeConverters[str] = converter
    end
    nothing
  end

  function specialprotoread(iob, instanceType::Type, Tstr::ASCIIString)
    it = instanceType{ ProtoTest.regTypes[Tstr] }()
    nothingspecialprotoread(iob, it)
  end

  function runOutsideFunction(f::Function)
    f()
    println("its done")
    nothing
  end

end

using MyMod, ProtoTest

mydatas = MyTemplType{MyType}(rand(2),1,"test01",MyType(rand(3),2, "testT"))

iob = protostring(mydatas)


registerType("MyType", MyType)
dd = specialprotoread(iob, MyTemplType, "MyType")

@show dd
@show mydatas

@test dd.x[1] == mydatas.x[1]



# dd = specialprotoread(iob, MyTemplType, MyType)

# dd = nothingspecialprotoread(iob, MyTemplType{MyType}())

# runOutsideFunction(MyMod.thing)










#
