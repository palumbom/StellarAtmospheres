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
function row_for_species(df::DataFrame, species::String)
    # find the appropriate row
    ind = findfirst(df.Species .== species)
    return convert(Array, df[ind, :][2:end-1])
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
        header = ["Z", "Species", "Weight", "First"]
        df = CSV.read(dir * "nist_ioniz.txt", header=header, delim="\t",
                      ignorerepeated=true, silencewarnings=true)
        col_to_float!(df, names(df)[3:end]...)
        df.Species = strip.(df.Species)
    else
        header = ["Z", "Species", "Weight", "First", "Second", "Third"]
        df = CSV.read(dir * "nist_ioniz.txt", header=header, delim=" ",
                      ignorerepeated=true, silencewarnings=true)
        col_to_float!(df, names(df)[3:end]...)
    end
    return df
end

const dfi = tabulate_ionization(nist=true)
