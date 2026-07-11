module IncrInfrApproxMinDegreeExt

using AMD
import IncrementalInference: _ccolamd, _ccolamd!

# elseif ordering == :ccolamd
#   cons = zeros(SuiteSparse_long, length(adjMat.colptr) - 1)
#   cons[findall(x -> x in constraints, permuteds)] .= 1
#   p = Ccolamd.ccolamd(adjMat, cons)
#   @warn "Ccolamd is experimental in IIF at this point in time."

const KNOBS = 20
const STATS = 20



function _ccolamd!(
  n_row, #SuiteSparse_long,
  A::AbstractVector{T}, # SuiteSparse_long},
  p::AbstractVector, # SuiteSparse_long},
  knobs::Union{Ptr{Nothing}, Vector{Float64}},
  stats::AbstractVector, #{SuiteSparse_long},
  cmember::Union{Ptr{Nothing}, <:AbstractVector}, #{SuiteSparse_long}},
) where T
  n_col = length(p) - 1

  if length(stats) != STATS
    error("stats must hcae length $STATS")
  end
  if isa(cmember, Vector) && length(cmember) != n_col
    error("cmember must have length $n_col")
  end

  Alen = AMD.ccolamd_l_recommended(length(A), n_row, n_col)
  resize!(A, Alen)

  for i in eachindex(A)
    A[i] -= 1
  end
  for i in eachindex(p)
    p[i] -= 1
  end
  # BSD-3 clause, (c) Davis, Rajamanickam, Larimore
  # https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/f98e0f5a69acb6a3fb19703ff266100d43491935/LICENSE.txt#L153
  err = AMD.ccolamd_l(
    n_row,
    n_col,
    Alen,
    A,
    p,
    knobs,
    stats,
    cmember
  )

  if err == 0
    AMD.ccolamd_l_report(stats)
    error("call to ccolamd return with error code $(stats[4])")
  end

  for i in eachindex(p)
    p[i] += 1
  end

  pop!(p) # remove last zero from pivoting vector
  return p
end

function _ccolamd!(
  n_row,
  A::AbstractVector{T1}, #SuiteSparse_long},
  p::AbstractVector{<:Real}, # {SuiteSparse_long},
  cmember::Union{Ptr{Nothing}, <:AbstractVector{T}}, # SuiteSparse_long
) where {T1<:Real, T<:Integer}
  n_col = length(p) - 1

  if length(cmember) != n_col
    error("cmember must have length $n_col")
  end

  Alen = AMD.ccolamd_l_recommended(length(A), n_row, n_col)
  resize!(A, Alen)
  stats = zeros(T1, STATS)
  return _ccolamd!(n_row, A, p, C_NULL, stats, cmember)
end

# function _ccolamd!(
#   n_row,
#   A::AbstractVector{T}, # ::Vector{SuiteSparse_long},
#   p::AbstractVector, # ::Vector{SuiteSparse_long},
#   constraints = zeros(T,length(p) - 1), # SuiteSparse_long, 
# ) where T
#   n_col = length(p) - 1
#   return _ccolamd!(n_row, A, p, constraints)
# end

_ccolamd(n_row,A,p,constraints) = _ccolamd!(n_row, copy(A), copy(p), constraints)
_ccolamd(biadjMat, constraints) = _ccolamd(size(biadjMat, 1), biadjMat.rowval, biadjMat.colptr, constraints)



end