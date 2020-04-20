# import special function needed for voigt
import SpecialFunctions.erfcx

using PyCall
apm = pyimport("astropy.modeling.functional_models")

"""
    voigt(u, a)

Return a voigt profile. Gray uses the Hjerting function (see
Gray Eq. 11.46 and text in vicinity). This expression (source??) is from
the complementary error function, which you can express the Faddeeva
function in terms of, which you can express the Voigt function in terms off
(woof). But this is much faster than manually computing and convolving
Gaussian and Lorentzian profiles (you'll just have to believe me).
"""
function voigt(u::T, a::T, ΔνD::T) where T<:AF
    return real(faddeeva(u + im * a))/(sqrt(π) * ΔνD)
end

function hjerting(u::T, a::T) where T<:AF
    return real(faddeeva(u + im * a))
end

function faddeeva(x::Complex{T}) where T<:AF
    return erfcx(-im * x)
end

"""
    c

See Gray Eq. 11.46
"""
function calc_α(λ::T, temp::T, Pe::T, Pg::T, ξ::T, line::LineParams) where T<:AF
    return line(λ, temp, Pe, Pg, ξ)
end

# function calc_α(λ::T, temp::AA{T,1}, Pe::AA{T,1}, Pg::AA{T,1}, ξ::AA{T,1}, line::LineParams) where T<:AF
#     map(x -> line(λ, x,))
#     return line(λ, temp, Pe, Pg, ξ)
# end

function calc_α(λ::AA{T,1}, temp::T, Pe::T, Pg::T, ξ::T, line::LineParams) where T<:AF
    return line.(λ, temp, Pe, Pg, ξ)
end

"""
    c

See Gray Eq. 11.46
"""
function (line::LineParams)(λ::T, temp::T, Pe::T, Pg::T, ξ::T) where T<:AF
    # we need oscillator strength & doppler factor
    f = calc_f(line)
    ΔνD = calc_ΔνD(temp, ξ, line)

    # first calculate u and a from params provided
    u = calc_u(λ, ΔνD, line)
    a = calc_a(temp, Pe, Pg, ΔνD, line)

    # put it all together and return
    factor = π * e^2/(me * c) * f
    return factor * voigt(u, a, ΔνD)
end

"""

"""
function calc_u(λ::T, ΔνD::T, line::LineParams) where T<:AF
    # wavelength must be in cm here
    return (λ2ν(λ*1e-8) - line.ν₀)/(ΔνD)
end

"""

"""
function calc_a(temp::T, Pe::T, Pg::T, ΔνD::T, line::LineParams) where T<:AF
    γ = calc_γ(temp, Pe, Pg, line)
    return γ / (4.0 * π * ΔνD)
end

"""
    calc_f()

Compute the oscillator strength as in Gray Eq. 11.12
"""
function calc_f(line::LineParams) where T<:AF
    return 1.884e-15 * line.λ₀^2 * (line.gu/line.gl) * line.A[1]
end

"""
    calc_γ()

Lorentzian widths add naively. See Gray Eq. 11.44 and subsequent text.
"""
function calc_γ(temp::T, Pe::T, Pg::T, line::LineParams) where T<:AF
    γn = calc_γn(line)
    γ4 = calc_γ4(temp, Pe, line)
    γ6 = calc_γ6(temp, Pg, line)
    return γn + γ4 + γ6
end

"""
    calc_γn()

Natural broadening damping factor. Gray Eq. 11.15 and subsequent text.
"""
function calc_γn(line::LineParams)
    return 4.0 * π * sum(line.A)
end

"""
    calc_γ4()

Stark broadening damping factor. Gray Eq. 11.27.
"""
function calc_γ4(temp::T, Pe::T, line::LineParams) where T<:AF
    return exp10(19.0 + (2.0/3.0) * line.logC4 + log10(Pe) - (5.0/6.0) * log10(temp))
end

"""
    calc_γ6()

Van der Waals broadening damping factor. Gray Eq. 11.29.
"""
function calc_γ6(temp::T, Pg::T, line::LineParams) where T<:AF
    return exp10(20.0 + 0.4 * line.logC6 + log10(Pg) - 0.7 * log10(temp))
end

# doppler wavelength
calc_ΔλD(temp::T, ξ::T, line::LineParams) where T<:AF = calc_ΔλD(line.λ₀, temp, line.m, ξ)
function calc_ΔλD(λ₀::T, temp::T, m::T, ξ::T) where T<:AF
    return (λ₀/c) * sqrt((2.0*kB*temp/m) + ξ^2.0)
end

calc_ΔνD(temp::T, ξ::T, line::LineParams) where T<:AF = calc_ΔνD(line.ν₀, temp, line.m, ξ)
function calc_ΔνD(ν₀::T, temp::T, m::T, ξ::T) where T<:AF
    return (ν₀/c) * sqrt((2.0*kB*temp/m) + ξ^2.0)
end
