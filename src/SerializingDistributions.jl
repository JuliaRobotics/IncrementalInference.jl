
# TODO stop-gap string storage of Distrubtion types, should be upgraded to more efficient storage
function normalfromstring(str::AS) where {AS <: AbstractString}
  meanstr = match(r"μ=[+-]?([0-9]*[.])?[0-9]+", str).match
  mean = split(meanstr, '=')[2]
  sigmastr = match(r"σ=[+-]?([0-9]*[.])?[0-9]+", str).match
  sigma = split(sigmastr, '=')[2]
  Normal{Float64}(parse(Float64,mean), parse(Float64,sigma))
end

function mvnormalfromstring(str::AS) where {AS <: AbstractString}
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

function categoricalfromstring(str::AS)::Distributions.Categorical where {AS <: AbstractString}
  # pstr = match(r"p=\[", str).match
  psubs = split(str, '=')[end]
  psubs = split(psubs, '[')[end]
  psubsub = split(psubs, ']')[1]
  pw = split(psubsub, ',')
  p = parse.(Float64, pw)
  return Categorical(p ./ sum(p))
end

function extractdistribution(str::AS)::Union{Void, Distributions.Distribution} where {AS <: AbstractString}
  # TODO improve use of multidispatch and packing of Distribution types
  if str == ""
    return nothing
  elseif (ismatch(r"Normal", str) && !ismatch(r"FullNormal", str))
    return normalfromstring(str)
  elseif ismatch(r"FullNormal", str)
    return mvnormalfromstring(str)
  elseif ismatch(r"Categorical", str)
    return categoricalfromstring(str)
  elseif ismatch(r"KDE:", str)
    return convert(KDE.BallTreeDensity, str)
  elseif ismatch(r"AliasingScalarSampler", str)
    return convert(AliasingScalarSampler, str)
  else
    error("Don't know how to extract distribution from str=$(str)")
  end
end
