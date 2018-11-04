module Ccolamd

using SparseArrays
using SuiteSparse.CHOLMOD: SuiteSparse_long

const KNOBS = 20
const STATS = 20

function recommended(nnz::SuiteSparse_long, n_row::SuiteSparse_long, n_col::SuiteSparse_long)
    Alen = ccall((:ccolamd_l_recommended, :libccolamd), Csize_t,
        (SuiteSparse_long, SuiteSparse_long, SuiteSparse_long),
        nnz, n_row, n_col)
    if Alen == 0
        error("")
    end
    return Alen
end

function set_defaults(knobs::Vector{Float64})
    if length(knobs) != KNOBS
        error("")
    end
    ccall((:ccolamd_set_defaults, :libccolamd), Nothing, (Ptr{Cdouble},), knobs)
    return knobs
end

function ccolamd!(n_row::SuiteSparse_long, A::Vector{SuiteSparse_long},
    p::Vector{SuiteSparse_long}, knobs::Union{Ptr{Nothing},Vector{Float64}}, stats::Vector{SuiteSparse_long},
    cmember::Union{Ptr{Nothing},Vector{SuiteSparse_long}})

    n_col = length(p) - 1

    if length(stats) != STATS
        error("stats must hcae length $STATS")
    end
    if isa(cmember, Vector) && length(cmember) != n_col
        error("cmember must have length $n_col")
    end

    Alen = recommended(length(A), n_row, n_col)
    resize!(A, Alen)

    for i in eachindex(A)
        A[i] -= 1
    end
    for i in eachindex(p)
        p[i] -= 1
    end
    err = ccall((:ccolamd_l, :libccolamd), SuiteSparse_long,
        (SuiteSparse_long, SuiteSparse_long, SuiteSparse_long, Ptr{SuiteSparse_long},
         Ptr{SuiteSparse_long}, Ptr{Cdouble}, Ptr{SuiteSparse_long}, Ptr{SuiteSparse_long}),
        n_row, n_col, Alen, A, p, knobs, stats, cmember)
    if err == 0
        report(stats)
        error("call to ccolamd return with error code $(stats[4])")
    end

    for i in eachindex(p)
        p[i] += 1
    end

    pop!(p) # remove last zero from pivoting vector
    return p
end

function ccolamd!(n_row, A::Vector{SuiteSparse_long}, p::Vector{SuiteSparse_long}, cmember::Union{Ptr{Nothing},Vector{SuiteSparse_long}})
    n_col = length(p) - 1

    if length(cmember) != n_col
        error("cmember must have length $n_col")
    end

    Alen = recommended(length(A), n_row, n_col)
    resize!(A, Alen)
    stats = zeros(SuiteSparse_long, STATS)
    return ccolamd!(n_row, A, p, C_NULL, stats, cmember)
end

function ccolamd!(n_row, A::Vector{SuiteSparse_long}, p::Vector{SuiteSparse_long})
    n_col = length(p) - 1
    cmember = zeros(SuiteSparse_long, length(p) - 1)
    return ccolamd!(n_row, A, p, cmember)
end

ccolamd(n_row, A::Vector{SuiteSparse_long}, p::Vector{SuiteSparse_long}) =
    ccolamd!(n_row, copy(A), copy(p))

ccolamd(A::SparseMatrixCSC) = ccolamd(size(A, 1), A.rowval, A.colptr)

report(stats::Vector{SuiteSparse_long}) =
    ccall((:ccolamd_l_report, :libccolamd), Nothing, (Ptr{SuiteSparse_long},), stats)

end #module
