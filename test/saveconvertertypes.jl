using Base: Test


module  Dependency
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
  using Dependency
  import Dependency: convert

  export T1, PackedT1, convertsave

  type T1 <: abst  end
  type PackedT1 <: pabst  end
end


using Extend

@test convertsave(T1()) == Extend.PackedT1








#
