# modules for reading/writing table data
using CSV
using DataFrames

datdir = (@__DIR__) * "/../data/"

# extend tryparse function
import Base.tryparse
function tryparse(Float64, x::Float64)
    return x
end

function tryparse(Float64, x::Missing)
    return NaN
end

function col_to_float!(df::DataFrame, colnames::T...) where T<:Symbol
    for col in colnames
        df[!, col] = tryparse.(Float64, df[!, col])
        df[isnothing.(df[!, col]), col] .= NaN
    end
    return df
end

"""

"""
function tabulate_partition(dir::String=datdir)
    @assert isdir(dir)

    # read in the files & # make sure numbers are floats
    df = CSV.read(dir * "partit.txt", header=0)
    col_to_float!(df, names(df)[2:end]...)
    return df
end

function tabulate_ionization(dir::String=datdir)
    @assert isdir(dir)

    # read in the files & # make sure numbers are floats
    df = CSV.read(dir * "ioniz.txt", header=0, delim=" ")
    col_to_float!(df, names(df)[2:end]...)
    return df
end
