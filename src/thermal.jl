# define planck functions
"""
    Bν(ν, T)

Compute the specific intensity at frequency ν of a blackbody with temperature T.
"""
Bν(ν::t, T::t) where t<:Real = (2.0*h*ν^3.0/c^2.0) * one(t)/(exp(h*ν/(kB*T)) - one(t))

"""
    Bλ(λ, T)

Compute the specific intensity at wavelength λ of a blackbody with temperature T.
"""
Bλ(λ::t, T::t) where t<:Real = (2.0*h*c^2.0/λ^5.0) * one(t)/(exp(h*c/(λ*kB*T)) - one(t))

"""
    Bν̃(ν̃, T)

Compute the specific intensity at wavenumber ν̃ of a blackbody with temperature T.
"""
Bν̃(ν̃::t, T::t) where t<:Real = (2.0*h*c^2.0*ν̃^3.0) * one(t)/(exp(h*c*ν̃/(kB*T)) - one(t))

SB(T::t) where t<:Real = σB * T^4

"""

"""
function Tτ(τ::T; Teff::T=NaN) where T<:AF
    @assert !isnan(Teff)
    return Teff * (0.75 * (τ + (2.0/3.0)))^0.25
end

"""

"""
function P_from_nkt(n::T, temp::T) where T<:AF
    return n * kB * temp
end

"""

"""
function ΦT(;temp::T=NaN, species::String="", ion::Symbol=:First) where T<:AF
    @assert !isnan(temp)
    @assert !isempty(species)
    return ΦT(temp, species, ion=ion)
end

function ΦT(temp::T, species::String; ion::Symbol=:First) where T<:AF
    if !(species in dfi.Species)
        println(species)
    end
    @assert species in dfp.Species

    # figure out what the next ionization is
    if species == "H"
        ionspec = "HII"
    elseif species == "H-"
        ionspec = "H"
    else
        ionspec = species * "+"
    end

    # now do the math
    θ = temp_to_theta(temp)
    U0 = exp10.(calc_partition(temp, species))
    U1 = exp10.(calc_partition(temp, ionspec))
    return 1.2020e9 * U1 * θ^(-5.0/2.0) * exp10(-θ * χ(species, ion)) / U0
end

"""

Calculate the neutral fraction of an element. Only consider 1st ionization.
"""
function neutral_fraction(temp::T, Pe::T, species::String) where T<:AF
    frac = Pe/ΦT(temp, species)
    return frac/(frac + one(T))
end

function species_fraction(temp::T, Pe::T, species::String; ion::String="Zeroth") where T<:AF
    @assert ion in ["Zeroth", "First", "Second"]
    if ((species == "H") | (species == "He"))
        @assert ion in ["Zeroth", "First"]
        n1n0 = ΦT(temp, "H") / Pe
        if ion == "Zeroth"
            return (one(T)/n1n0)/((one(T)/n1n0) + one(T))
        elseif ion == "First"
            return one(T)/((one(T)/n1n0) + one(T))
        end
    end

    # general case
    n1n0 =  ΦT(temp, species) / Pe
    n2n1 =  ΦT(temp, species * "+") / Pe

    if ion == "Zeroth"
        return (one(T)/n1n0)/((one(T)/n1n0) + one(T) + n2n1)
    elseif ion == "First"
        return one(T)/((one(T)/n1n0) + one(T) + n2n1)
    elseif ion == "Second"
        return n2n1/((one(T)/n1n0) + one(T) + n2n1)
    end
end

"""

"""
function calc_SE(λ::T, temp::T) where T<:AF
    θ = temp_to_theta(temp)
    χλ = 1.2398e4/λ
    return one(T) - exp10(-χλ * θ)
end
