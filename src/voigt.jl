import SpecialFunctions.erfcx

"""

"""
function voigt(u::T, a::T) where T<:AF
    x = (u + im * a)
    return real(erfcx(-im * x))
end

# functions for u and a
calc_Δλ(v_R::T, λ::T) where T<:AF = v_R * λ / c
calc_Δν(v_R::T, ν::T) where T<:AF = v_R * ν / c

function calc_ΔλD(λ₀::T, temp::T, μ::T) where T<:AF
    return 4.301e-7 * λ₀ * sqrt(temp/μ)
end

function calc_ΔνD(ν₀::T, temp::T, μ::T) where T<:AF
    return 4.301e-7 * ν₀ * sqrt(temp/μ)
end

function calc_ΔλD(λ₀::T, temp::T, m::T, ξ::T) where T<:AF
    return (λ₀/c) * sqrt((2.0*kB*temp/m) + ξ^2.0)
end

function calc_ΔνD(ν₀::T, temp::T, m::T, ξ::T) where T<:AF
    return (ν₀/c) * sqrt((2.0*kB*temp/m) + ξ^2.0)
end
