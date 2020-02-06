BG = Base.Generator
SRL = StepRangeLen
BTP = Base.TwicePrecision

"""
    logspace(a, b; length=NaN)

Compute logarithmically-spaced bins between a and b. Returns a generator.
"""
function logspace(a::T, b::T; length::Int=NaN) where T<:Real
    @assert b > a
    @assert !isnan(length)
    # return exp10.(range(a, stop=b, length=len))
    return (exp10(x) for x in range(log10(a), log10(b), length=length))
end

"""
    getindex(a, n)

Return the nth index of generator a, specifically for generator of type
returned by logspace
"""
function getindex(a::BG{SRL{T,BTP{T},BTP{T}},typeof(exp10)}, n::Int) where T<:Real
    return first(Iterators.drop(a, n-1))
end
