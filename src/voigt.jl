import SpecialFunctions.erfcx

"""

"""
function voigt(u::T, a::T) where T<:AF
    x = (u + im * a)
    return real(erfcx(-im * x))
end

"""

"""
function calc_u(x::T, λ₀::T, temp::T, m::T, ξ::T) where T<:AF
    return (x - λ₀)/calc_ΔλD(λ₀, temp, m, ξ)
end

"""

"""
function calc_a(Pe::T, temp::T, line::LineParams) where T<:AF
    γn = calc_γn()
    γ4 = calc_γ4()
    γ6 = calc_γ6()
    return
end

# natural broadening
function calc_γn()

    return
end

# stark broadening
function calc_γ4(Pe::T, temp::T, line::LineParams) where T<:AF
    return 19.0 + (2.0/3.0) * line.logC4 + log10(Pe) - (5.0/6.0) * log10(temp)
end

# van der waals broadening
function calc_γ6(Pg::T, temp::T, line::LineParams)
    return 20.0 + 0.4 * line.logC6 + log10(Pg) - 0.7 * log10(temp)
end

# doppler wavelength
function calc_ΔλD(λ₀::T, temp::T, m::T, ξ::T) where T<:AF
    return (λ₀/c) * sqrt((2.0*kB*temp/m) + ξ^2.0)
end

function calc_ΔνD(ν₀::T, temp::T, m::T, ξ::T) where T<:AF
    return (ν₀/c) * sqrt((2.0*kB*temp/m) + ξ^2.0)
end
