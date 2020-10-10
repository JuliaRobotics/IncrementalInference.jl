
# TODO stop-gap string storage of Distrubtion types, should be upgraded to more efficient storage
function normalfromstring(str::AbstractString)
  meanstr = match(r"μ=[+-]?([0-9]*[.])?[0-9]+", str).match
  mean = split(meanstr, '=')[2]
  sigmastr = match(r"σ=[+-]?([0-9]*[.])?[0-9]+", str).match
  sigma = split(sigmastr, '=')[2]
  Normal{Float64}(parse(Float64,mean), parse(Float64,sigma))
end

function mvnormalfromstring(str::AbstractString)
  means = split(split(split(str, 'μ')[2],']')[1],'[')[end]
  mean = Float64[]
  for ms in split(means, ',')
    push!(mean, parse(Float64, ms))
  end
  sigs = split(split(split(str, 'Σ')[2],']')[1],'[')[end]
  sig = Float64[]
  for ms in split(sigs, ';')
    for m in split(ms, ' ')
      length(m) > 0 ? push!(sig, parse(Float64, m)) : nothing
    end
  end
  len = length(mean)
  sigm = reshape(sig, len,len)
  MvNormal(mean, sigm)
end

function categoricalfromstring(str::AbstractString)
  # pstr = match(r"p=\[", str).match
  psubs = split(str, '=')[end]
  psubs = split(psubs, '[')[end]
  psubsub = split(psubs, ']')[1]
  pw = split(psubsub, ',')
  p = parse.(Float64, pw)
  return Categorical(p ./ sum(p))
end



function extractdistributionJson(packedSamplable::Dict{Symbol,String})
  # TODO improve use of multidispatch and packing of Distribution types
  @show packedSamplable

  error("Not implemented yet")
end


# NOTE part of new effort to overhaul the SamplableBelief serialization approach
# maybe it all becomes a JSON struct sort of thing in the long run.
StringThemSamplableBeliefs = Union{Normal, MvNormal, Categorical, DiscreteNonParametric, BallTreeDensity, AliasingScalarSampler}
convert(::Type{<:PackedSamplableBelief}, obj::StringThemSamplableBeliefs) = string(obj)


function convert(::Type{<:SamplableBelief}, str::Union{<:PackedSamplableBelief,<:AbstractString})
  # extractdistribution(str::AS) where {AS <: AbstractString}
  # TODO improve use of multidispatch and packing of Distribution types
  # TODO use startswith
  if str == ""
    return nothing
  elseif occursin(r"SamplableTypeJSON", str)
    # TODO this is the new direction for serializing (pack/unpack) of <:Samplable objects
    # NOTE uses intermediate consolidation keyword search pattern `SamplableTypeJSON`
    props = JSON2.read(str, Dict{Symbol,String})
    return extractdistributionJson(props)
  elseif startswith(str, "DiagNormal")
    # Diags are internally squared, so only option here is to sqrt on input.
    return mvnormalfromstring(str)
  elseif (occursin(r"Normal", str) && !occursin(r"FullNormal", str))
    return normalfromstring(str)
  elseif occursin(r"FullNormal", str)
    return mvnormalfromstring(str)
  elseif occursin(r"Categorical", str)
    return categoricalfromstring(str)
  elseif occursin(r"DiscreteNonParametric", str)
    return categoricalfromstring(str)
  elseif occursin(r"KDE:", str)
    return convert(BallTreeDensity, str)
  elseif occursin(r"AliasingScalarSampler", str)
    return convert(AliasingScalarSampler, str)
  else
    error("Don't know how to extract distribution from str=$(str)")
  end
end


function convert(::Type{<:PackedSamplableBelief}, obj::SamplableBelief)
  # FIXME must use string, because unpacking templated e.g. PackedType{T} has problems, see DFG #668
  string(obj)
end



#