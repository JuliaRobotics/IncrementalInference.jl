# test module wrapping

using Base.Test

module First


export solve, InferenceType, Container

abstract InferenceType

type Container
  col::Dict{Int, Function}
end

function registerCallback!(col::Container, id::Int, f::Function)
  col.col[id] = f
end

function solve(ctl::Container, id::Int, val::InferenceType)
  if false
    # doesnt work
    evalPotential(val)
  else
    evalPotential = ctl.col[id]
    evalPotential(val)
  end
end

end




module Second
importall First
export SecondType, SecondAgain, evalPotential, solve, registerCallback!, Container

type SecondType <: First.InferenceType
  x::Int
end
type SecondAgain <: First.InferenceType
  x::Int
end


function evalPotential(x::SecondType)
  println("evalPotential sees $(x)")
  return x.x
end

function evalPotential(x::SecondAgain)
  println("evalPotential also sees $(x)")
  return x.x
end


end


# Okay lets see
using Second

Col = Container(Dict{Int, Function}())

First.registerCallback!(Col, 1, Second.evalPotential)

stA = 1.0
saA = 3.0

st = SecondType(stA)
sa = SecondAgain(saA)


@test solve(Col, 1, st) == stA
@test solve(Col, 1, sa) == saA

@elapsed evalPotential(st)
t1 = @elapsed evalPotential(st)
t2 = @elapsed solve(Col, 1, sa)

println("Check the speed is reasonable")
@test t2 <15.0*t1 # should actually be about equal, slack for threading uncertainty





#
