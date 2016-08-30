
type MyType
   x::Array{Float64,1}
   str::ASCIIString
   MyType() = new()
   MyType(a::Array{Float64,1}, b::ASCIIString) = new(a,b)
end

mkexp(str,x...) = :($(parse(str))($(x...)))
# mkexp(ty,str,x...) = :(ty$(parse(str))($(x...)))
strT = "MyType"
@show ex1 = mkexp(strT);
@show typeof(ex1);
mt1 = eval(ex1)

ex2 = mkexp(strT);
mt1 = eval(ex2)
