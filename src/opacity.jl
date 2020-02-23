# modules for reading/writing table data
using CSV
using Dierckx
using DataFrames

# directory path for the data
datdir = (@__DIR__) * "/../data/"

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

"""

"""
function tabulate_ionization(dir::String=datdir)
    @assert isdir(dir)

    # read in the files & # make sure numbers are floats
    header = ["A", "Element", "Weight", "First", "Second", "Third"]
    df = CSV.read(dir * "ioniz.txt", header=header, delim=" ",
                  ignorerepeated=true, silencewarnings=true)
    col_to_float!(df, names(df)[3:end]...)
    return df
end

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

"""
function calc_partition(temp::T, species::String) where T<:Real
    df = tabulate_partition()

    @assert species in df.Species
    @assert temp >= 0

    # check for simple species
    if species == "H-"
        return 1.0
    elseif species == "HII"
        return 1.0
    end

    # theta data
    θ = range(0.2, 2.0, step=0.2)
    P = row_for_species(df, species)

    # now do the interpolation
    spl = Spline1D(θ, P)
    return spl(temp_to_theta(temp))
end
