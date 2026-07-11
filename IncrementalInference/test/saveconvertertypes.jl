using Test


module Dependency
  import Base: convert
  export abst, pabst, convert, convertsave

  abstract type abst end
  abstract type pabst end

  convert(::Type{P}, ::T) where {P <: pabst, T <: abst} =
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


@testset "Ensure converter types can be run from extending namespaces..." begin

@test convertsave(T1()) == Extend.PackedT1

end



#
