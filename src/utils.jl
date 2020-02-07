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
    return range(log(a), log(b), length=length)
end
