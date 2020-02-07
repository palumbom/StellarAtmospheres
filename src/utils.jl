BG = Base.Generator
SRL = StepRangeLen
BTP = Base.TwicePrecision

"""
    logspace(a, b; length=NaN)

Compute logarithmically-spaced bins between a and b.
"""
function logspace(a::T, b::T; length::Int=NaN) where T<:Real
    @assert b > a
    @assert !isnan(length)
    return exp10(range(log10(a), log10(b), length=length))
end

asum(a::AbstractVector) = asum(a, dims=1)
function asum(a::AbstractArray{T,N}; dims::Integer) where {T,N}
    1 <= dims <= N || throw(ArgumentError("dimension $dims out of range (1:$N)"))

    r = axes(a)
    r0 = ntuple(i -> i == dims ? UnitRange(1, last(r[i]) - 1) : UnitRange(r[i]), N)
    r1 = ntuple(i -> i == dims ? UnitRange(2, last(r[i])) : UnitRange(r[i]), N)

    return view(a, r1...) .+ view(a, r0...)
end
