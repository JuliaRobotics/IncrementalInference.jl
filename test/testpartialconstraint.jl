# partial constraint development

using IncrementalInference, Distributions


type DevelopPartial <: IncrementalInference.FunctorSingleton
  x::Distribution
  partial::Tuple
end
getSample(dpl::DevelopPartial, N::Int=1) = (rand(dpl.x, N), )


dp = DevelopPartial(Normal(),(1,))



















#
