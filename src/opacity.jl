# modules for reading/writing table data
using CSV
using Dierckx
using DataFrames

# directory path for the data
datdir = abspath((@__DIR__) * "/../data/")

# extend tryparse function
import Base.tryparse
function tryparse(Float64, x::Float64)
    return x
end

function tryparse(Float64, x::Missing)
    return NaN
end

function tryparse(Float64, x::Nothing)
    return NaN
end

"""
    col_to_float!(df, colnames)

Replace strings, missing, and nothing values in columns that should be floats.
"""
function col_to_float!(df::DataFrame, colnames::T...) where T<:Symbol
    for col in colnames
        df[!, col] = tryparse.(Float64, df[!, col])
        df[isnothing.(df[!, col]), col] .= NaN
        df[ismissing.(df[!, col]), col] .= NaN
        df[!, col] = convert.(Float64, df[!,col])
        df[!, col] = replace(df[!,col], -0.0=>NaN)
    end
    return df
end

"""

"""
function tabulate_partition(dir::String=datdir)
    @assert isdir(dir)

    # read in the files & # make sure numbers are floats
    num = string.(range(0.2, 2.0, step=0.2))
    header = ["Species", num..., "logg0"]
    df = CSV.read(dir * "partit.txt", header=header, silencewarnings=true)
    col_to_float!(df, names(df)[2:end]...)
    return df
end

const dfp = tabulate_partition()

"""

"""
function tabulate_ionization(dir::String=datdir; nist::Bool=true)
    @assert isdir(dir)

    # read in the files & # make sure numbers are floats
    if nist
        header = ["A", "Species", "Weight", "First"]
        df = CSV.read(dir * "nist_ioniz.txt", header=header, delim="\t",
                      ignorerepeated=true, silencewarnings=true)
        col_to_float!(df, names(df)[3:end]...)
        df.Species = strip.(df.Species)
    else
        header = ["A", "Species", "Weight", "First", "Second", "Third"]
        df = CSV.read(dir * "nist_ioniz.txt", header=header, delim=" ",
                      ignorerepeated=true, silencewarnings=true)
        col_to_float!(df, names(df)[3:end]...)
    end
    return df
end

const dfi = tabulate_ionization(nist=true)

"""

"""
function temp_to_theta(temp::T) where T<:Real
    return 5040.0/(temp)
end

function theta_to_temp(θ::T) where T<:Real
    return 5040.0/θ
end

"""

"""
function row_for_species(df::DataFrame, species::String)
    # find the appropriate row
    ind = findfirst(df.Species .== species)
    return convert(Array, df[ind, :][2:end-1])
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
function χ(n::Int, species::String; ion::Symbol=:First)
    return χ(species, ion) - ((h*R*c)/n^2) * eV
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
    ϴ = temp_to_theta(temp)
    λR = λc* R * 1e-8 # convert Å to cm
    return one(T) + (0.3456/λR^(1.0/3.0)) * (loge/(θ * 1.2398e4/λ) + 0.5)
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
        sm += g_bf(λ, nn) * exp10(-ϴ * χ1) / nn^3
    end

    # second term on opacity
    t2 = (loge/(2.0 * θ * χ1)) * (exp10(-χ3 * θ) - exp10(-θ * χ1))
    return α0 * λ^(3.0) * (sm + t2)
end
