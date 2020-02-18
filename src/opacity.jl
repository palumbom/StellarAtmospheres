# modules for reading/writing table data
using CSV
using DataFrames

datdir = (@__DIR__) * "/../data/"

# extend tryparse function
import Base.tryparse
function tryparse(Float64, x::Float64)
    return x
end

"""

"""
function tabulate_partition(dir::String=datdir)
    @assert isdir(dir)

    # read in the files
    df = CSV.read(dir * "partit.txt", header=0)

    # loop over rows
    for i in 2:length(names(df))
        df[!, names(df)[i]] = tryparse.(Float64, df[!, names(df)[i]])
    end
    return df
end
