# define composite type to hold line characteristics
struct LineParams{T<:AF}
    λ₀::T       # rest frame line center
    ν₀::T       # line center but ~frequency~
    A::AA{T,1}  # Einstein A coefficients
    m::T        # mass of atomic in cgs
    gu::T       # statistical weight of upper level
    gl::T       # statistical weight of lower level
    logC4::T    # stark interaction constant
    logC6::T    # van der waals interaction constant
end

# make an outer constructor so you can't mess up arg order
function LineParams(;λ₀=NaN, A=[NaN], m=NaN, gu=NaN, gl=NaN, logC4=NaN, logC6=NaN)
    ν₀ = λ2ν(λ₀)
    return LineParams(λ₀, ν₀, A, gu, gl, logC4, logC6)
end

"""
    calc_f()

Compute the oscillator strength as in Gray Eq. 11.12
"""
function calc_f(λ::T, line::LineParams) where T<:AF
    return 1.884e-15 * λ^2 * line.gu * line.A[1] / line.gl
end
