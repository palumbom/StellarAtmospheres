# import special function needed for voigt
import SpecialFunctions.erfcx

"""
    voigt(u, a)

Return a voigt profile. Gray calls this the Hjerting function (see
Gray Eq. 11.46 and text in vicinity. This expression (source??) is from
the complementary error function, which you can express the Faddeeva
function in terms of, which you can express the Voigt function in terms off
(woof). But this is much faster than manually computing and convolving
Gaussian and Lorentzian profiles (you'll just have to believe me).
"""
function voigt(u::T, a::T) where T<:AF
    x = (u + im * a)
    return real(erfcx(-im * x))
end

"""
    calc_α(x, Pe, Pg, temp, ξ, line)

See Gray Eq. 11.46
"""
function (line::LineParams)(x::T, Pe::T, Pg::T, temp::T, ξ::T) where T<:AF
    # first calculate u and a from params provided
    u = calc_u(x, temp, ξ, line)
    a = calc_a(Pe, Pg, temp, ξ, line)

    # now we need oscillator strength & doppler factor
    f = calc_f(x, line)
    ΔλD = calc_ΔλD(temp, ξ, line)

    # put it all together and return
    factor = (sqrt(π)* e^2/(line.m * c^2)) * (line.λ₀^2 * f/ΔλD)
    return factor * voigt(u, a)
end


"""

"""
function calc_u(x::T, temp::T, ξ::T, line::LineParams) where T<:AF
    return (x - line.λ₀)/calc_ΔλD(temp, ξ, line)
end

"""

"""
function calc_a(Pe::T, Pg::T, temp::T, ξ::T, line::LineParams) where T<:AF
    γ = calc_γ(Pe, Pg, temp, line)
    return (γ/(4π)) * (line.λ₀^2/c) * one(T)/calc_ΔλD(temp, ξ, line)
end

"""
    calc_γ()

Lorentzian widths add naively. See Gray Eq. 11.44 and subsequent text.
"""
function calc_γ(Pe::T, Pg::T, temp::T, line::LineParams) where T<:AF
    γn = calc_γn(line)
    γ4 = calc_γ4(Pe, temp, line)
    γ6 = calc_γ6(Pg, temp, line)
    return γn + γ4 + γ6
end

"""
    calc_γn()

Natural broadening damping factor. Gray Eq. 11.15 and subsequent text.
"""
function calc_γn(line::LineParams)
    # γu =
    return sum(line.A)
end

"""
    calc_γ4()

Stark broadening damping factor. Gray Eq. 11.27.
"""
function calc_γ4(Pe::T, temp::T, line::LineParams) where T<:AF
    return exp10(19.0 + (2.0/3.0) * line.logC4 + log10(Pe) - (5.0/6.0) * log10(temp))
end

"""
    calc_γ6()

Van der Waals broadening damping factor. Gray Eq. 11.29.
"""
function calc_γ6(Pg::T, temp::T, line::LineParams) where T<:AF
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
