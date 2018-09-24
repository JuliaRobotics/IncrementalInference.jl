using Base: Test


module Dependency
  using Compat
  import Base: convert
  export abst, pabst, convert, convertsave

  @compat abstract type abst end
  @compat abstract type pabst end

  convert{P <: pabst, T <: abst}(::Type{P}, ::T) =
          getfield(T.name.module, Symbol("Packed$(T.name.name)"))
  convertsave(t) = convert(pabst, t)
end

module Extend
  using Main.Dependency
  import Main.Dependency: convert

  export T1, PackedT1, convertsave

  mutable struct T1 <: abst  end
  mutable struct PackedT1 <: pabst  end
end


using Main.Extend

@test convertsave(T1()) == Extend.PackedT1








#
