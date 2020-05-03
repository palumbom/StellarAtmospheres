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
    return exp10.(range(log10(a), log10(b), length=length))
end

"""

Compute sum of adjacent array elements along specified dim.
"""
asum(a::AbstractVector) = asum(a, dims=1)
function asum(a::AbstractArray{T,N}; dims::Integer) where {T,N}
    1 <= dims <= N || throw(ArgumentError("dimension $dims out of range (1:$N)"))

    r = axes(a)
    r0 = ntuple(i -> i == dims ? UnitRange(1, last(r[i]) - 1) : UnitRange(r[i]), N)
    r1 = ntuple(i -> i == dims ? UnitRange(2, last(r[i])) : UnitRange(r[i]), N)
    return view(a, r1...) .+ view(a, r0...)
end

"""

Compute average between adjacent elements along specified dim.
"""
elav(a::AbstractVector) = elav(a, dims=1)
function elav(a::AbstractArray{T,N}; dims::Integer) where {T,N}
    1 <= dims <= N || throw(ArgumentError("dimension $dims out of range (1:$N)"))

    r = axes(a)
    r0 = ntuple(i -> i == dims ? UnitRange(1, last(r[i]) - 1) : UnitRange(r[i]), N)
    r1 = ntuple(i -> i == dims ? UnitRange(2, last(r[i])) : UnitRange(r[i]), N)
    return (view(a, r1...) .+ view(a, r0...)) ./ 2.0
end

"""Credit: traktofon @ https://discourse.julialang.org/t/findnearest-function/4143/4"""
function searchsortednearest(a::AbstractVector, x::Real) #where T<:Real
   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end

function searchsortednearest(x::Real, a::AbstractVector)
   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end
