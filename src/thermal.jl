# define planck functions
"""
    Bν(ν, T)

Compute the specific intensity at frequency ν of a blackbody with temperature T.
"""
Bν(ν::t, T::t) where t<:Real = (2.0*h*ν^3/c^2) * one(t)/(exp(h*ν/(kB*T)) - one(t))

"""
    Bλ(λ, T)

Compute the specific intensity at wavelength λ of a blackbody with temperature T.
"""
Bλ(λ::t, T::t) where t<:Real = (2.0*h*c^2/λ^5) * one(t)/(exp(h*c/(λ*kB*T)) - one(t))

"""
    Bν̃(ν̃, T)

Compute the specific intensity at wavenumber ν̃ of a blackbody with temperature T.
"""
Bν̃(ν̃::t, T::t) where t<:Real = (2.0*h*c^2*ν̃^3) * one(t)/(exp(h*c*ν̃/(kB*T)) - one(t))

SB(T::t) where t<:Real = σB * T^4


function Tτ(τ::T; Teff::T=NaN) where T<:Real
    @assert !isnan(Teff)
    return Teff * (0.75 * (τ + (2.0/3.0)))^0.25
end


# wien displacement law (cm and Hz and K)
λmax(T::t) where t<:Real = 0.290/T
νmax(T::t) where t<:Real = 5.88e10 * T

# convert λ to ν and ν to λ
λ2ν(λ::T) where T<:Real = c/λ
ν2λ(ν::T) where T<:Real = c/ν
λ2ν̃(λ::T) where T<:Real = one(T)/λ
ν̃2λ(ν̃::T) where T<:Real = one(T)/ν̃
ν̃2ν(ν̃::T) where T<:Real = c*ν̃
ν2ν̃(ν::T) where T<:Real = ν/c
