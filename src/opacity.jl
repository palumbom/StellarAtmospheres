"""

"""
function temp_to_theta(temp::T) where T<:Real
    return 5040.0/(temp)
end

function theta_to_temp(θ::T) where T<:Real
    return 5040.0/θ
end

"""
Calculate partition functions from Gray using his Table D.2 to interpolate
"""
function calc_partition(temp::T, species::String) where T<:AF
    @assert temp >= zero(T)

    # check for simple species
    if species == "H-"
        return one(T)
    elseif species == "HII"
        return one(T)
    end

    # now check for species in table
    @assert species in dfp.Species

    # theta data
    θs = range(0.2, 2.0, step=0.2)
    Ps = row_for_species(dfp, species)

    # now do the interpolation
    spl = Spline1D(θs, Ps)
    return spl(temp_to_theta(temp))
end

"""

"""
function χ(species::String, ion::Symbol)
    @assert ion in names(dfi)

    # immediately return for H- in eV
    if species == "H-"
        return 0.754195
    end

    # else read from the table
    return dfi[findfirst(dfi.Species .== species), ion]
end


"""

"""
function ΦT(;temp::T=NaN, species::String="", ion::Symbol=:First) where T<:AF
    @assert !isnan(temp)
    @assert !isempty(species)
    return ΦT(temp, species, ion=ion)
end

function ΦT(temp::T, species::String; ion::Symbol=:First) where T<:AF
    @assert species in dfi.Species
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
    U0 = exp10(calc_partition(temp, species))
    U1 = exp10(calc_partition(temp, ionspec))
    return 1.2020e9 * U1 * θ^(-5.0/2.0) * exp10(-θ * χ(species, ion)) / U0
end

"""

Gray Eq. 8.3
"""
function χ(n::Int, species::String)
    return χ(species, :First) - ((h*R*c)/n^2) * eV
end

"""

Gray Eq. 8.4. Returns value in cm^2 per neutral H atom
"""
function α_bf_H(λ::T, n::Int) where T<:AF
    return α0 * g_bf(λ, n) * (λ^3.0/n^5.0)
end

"""

Gray Eq. 8.5. Return value is dimensionless.
"""
function g_bf(λ::T, n::Int) where T<:AF
    λR = λ * R * 1e-8 # convert Å to cm
    return one(T) - (0.3456/λR^(1.0/3.0)) * (λR/n^2 - 0.5)
end

"""

Gray Eq. 8.6. Return value is dimensionless.
"""
function g_ff(λ::T, temp::T, Pe::T) where T<:AF
    θ = temp_to_theta(temp)
    λR = λ * R * 1e-8 # convert Å to cm
    return one(T) + ((0.3456/λR^(1.0/3.0)) * (loge/(θ * 1.2398e4/λ) + 0.5))
end

"""

Gray Eq. 8.8.
"""
function κ_H_bf(λ::T, temp::T, Pe::T) where T<:AF
    # tabulate some needed numbers
    λcm = λ * 1e-8
    θ = temp_to_theta(temp)
    n0 = floor(Int, sqrt(13.6/(eV*h*c/(λcm))) + one(T))  # dimensionless
    χ1 = χ("H", :First)
    χ3 = χ1 * (one(T) - (one(T)/(n0 + 3)^2))

    # now do some summing
    nf = n0 + 2
    nn = n0
    sm = 0.0

    while nn <= nf
        sm += g_bf(λ, nn) * exp10(-θ * χ(nn, "H")) / nn^3
        nn += 1
    end

    # second term on opacity
    t2 = (loge/(2.0 * θ * χ1)) * (exp10(-χ3 * θ) - exp10(-θ * χ1))
    return α0 * λ^(3.0) * (sm + t2)
end

"""

Gray Eq. 8.9.
"""
function α_ff_H(λ::T, temp::T, Pe::T) where T<:AF
    λcm = λ * 1e-8
    ν = λ2ν(λcm)
    num = 2.0 * h^(2.0) * e^(2.0) * R * ((2.0*mp)/(kB*temp))^(0.5)
    den = 3.0^(3.0/2.0) * π * mp^(3.0) * ν^3.0
    return num/den
end

"""

Gray Eq. 8.10.
"""
function κ_H_ff(λ::T, temp::T, Pe::T) where T<:AF
    θ = temp_to_theta(temp)
    χ1 = χ("H", :First)
    λcm = λ * 1e-8
    return α0 * λ^(3.0) * g_ff(λ, temp, Pe) * (loge/(2.0*θ*χ1)) * exp10(-θ*χ1)
end

"""

Gray Eq. 8.11. Units of 1e-18 cm^2/H- ion
"""
function α_bf_minus(λ::T) where T<:AF
    # short circuit return if long λ
    if λ > 16270.0
        return zero(T)
    end
    as = [1.99654, -1.18267e-5, 2.64243e-6, -4.40524e-10,
          3.23992e-14, -1.39568e-18, 2.78701e-23]
    return 1e-18 * sum([as[n] * λ^(n-1) for n in 1:length(as)])
end

"""

Gray Eq. 8.12.
"""
function κ_H_minus_bf(λ::T, temp::T, Pe::T) where T<:AF
    θ = temp_to_theta(temp)
    return 4.158e-10 * α_bf_minus(λ) * Pe * θ^(5.0/2.0) * exp10(0.754 * θ)
end

"""

Gray Eq. 8.13.
"""
function κ_H_minus_ff(λ::T, temp::T, Pe::T) where T<:AF
    θ = temp_to_theta(temp)
    logλ = log10(λ)

    # polynomial fit from Bell and Berrington
    f0 = -2.2763 - 1.6850*logλ + 0.76661*logλ^2.0 - 0.053346*logλ^3.0
    f1 = 15.2827 - 9.2846*logλ + 1.99381*logλ^2.0 - 0.142631*logλ^3.0
    f2 = -197.789 + 190.266*logλ - 67.9775*logλ^2.0 + 10.6913*logλ^3.0 - 0.625151*logλ^4.0
    return 1e-26 * Pe * exp10(f0 + f1*log10(θ) + f2*log10(θ)^2.0)
end

"""

Gray Eq. 8.18 abridged
"""
function κ_tot(λ::T, temp::T, Pe::T) where T<:AF
    # calculate each opacity
    κHbf = κ_H_bf(λ, temp, Pe)
    κHff = κ_H_ff(λ, temp, Pe)
    κHmbf = κ_H_minus_bf(λ, temp, Pe)
    κHmff = κ_H_minus_ff(λ, temp, Pe)

    # other quantities
    θ = temp_to_theta(temp)
    Φs = ΦT(temp, "H", ion=:First)
    χλ = 1.2398e4/λ
    return ((κHbf + κHff + κHmbf)*(one(T) - exp10(-χλ*θ)) + κHmff)/(one(T) + Φs/Pe)
end
