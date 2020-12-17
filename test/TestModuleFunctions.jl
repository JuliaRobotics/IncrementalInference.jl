# test module wrapping


# addprocs(1) # TODO reinsert once Julia 0.7 tests pass

using Test
using IncrementalInference

module First

export solve, AnotherInferenceType, Container

abstract type AnotherInferenceType end

mutable struct Container
  col::Dict{Symbol, Function}
  Container(::Nothing) = new()
  Container(;col=Dict{Symbol, Function}()) = new(col)
end

function registerCallback!(col::Container, fnc::Function)
  # get module
  # Symbol(string(m))
  m = Symbol(typeof(fnc).name.module)
  col.col[m] = fnc
end

function solve(ctl::Container, val::AnotherInferenceType) # m::Symbol,
  m = Symbol(typeof(val).name.module)
  if false
    # doesnt work
    evalPotential(val)
  else
    evalPotential = ctl.col[m]
    evalPotential(val)
  end
end

end


module SecondModule

using Main.First
export SecondType, SecondAgain, evalPotential, solve, registerCallback!, Container

struct SecondType <: Main.First.AnotherInferenceType
  x::Int
end
struct SecondAgain <: Main.First.AnotherInferenceType
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


# get module
# typeof(f).name.mt.name

# Okay lets see
using Main.SecondModule


@testset "out of module evalPotential..." begin


Col = First.Container()


# running on first process only
First.registerCallback!(Col, SecondModule.evalPotential)

CCol = deepcopy(Col)

stA = 1.0
saA = 3.0

st = SecondModule.SecondType(stA)
sa = SecondModule.SecondAgain(saA)

# Symbol(typeof(st).name.module)
@test First.solve(CCol, st) == stA
@test First.solve(CCol, sa) == saA


@elapsed SecondModule.evalPotential(st)
t1 = @elapsed SecondModule.evalPotential(st)
t2 = @elapsed solve(Col, sa)

println("Check the speed is reasonable")
@test t2 < 15.0*t1 # should actually be about equal, slack for threading uncertainty


end

# expand tests to include multiprocessor
# println("Tesing call of function on separate process...")
# fr = remotecall(evalPotential,procs()[2], st)
# @test fetch(fr) == stA
#
# fr = remotecall(solve,procs()[2], CCol,sa)
# @test fetch(fr) == saA
#
#
#
# println("Stopping all but first process...")
# rmprocs(procs()[2:end])
# @show procs()

#
