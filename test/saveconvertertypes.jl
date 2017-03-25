using Base: Test


module  Dependency
  import Base: convert
  export abst, pabst, convert, convertsave

  abstract abst
  abstract pabst

  convert{P <: pabst, T <: abst}(::Type{P}, ::T) =
          getfield(T.name.module, Symbol("Packed$(T.name.name)"))
  convertsave(t) = convert(pabst, t)
end

module Extend
  using Dependency
  import Dependency: convert

  export T1, PackedT1, convertsave

  type T1 <: abst  end
  type PackedT1 <: pabst  end

  convert{P <: pabst, T <: abst}(::Type{P}, ::T) =
        getfield(T.name.module, Symbol("Packed$(T.name.name)"))
end


using Extend

@test convertsave(T1()) == Extend.PackedT1








#
