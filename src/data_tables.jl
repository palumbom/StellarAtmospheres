# modules for reading/writing table data
using CSV
using Dierckx
using DataFrames
# using Interpolations
# using PyCall
# interp1d = pyimport("scipy.interpolate").interp1d

# directory path for the data
datdir = abspath((@__DIR__) * "/../data/")

# extend tryparse function
import Base.tryparse
function tryparse(Float64, x::Float64)
    return x
end

function tryparse(Float64, x::Int64)
    return float(x)
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
        df.Species = string.(strip.(df.Species))
    else
        header = ["Z", "Species", "Weight", "First", "Second", "Third"]
        df = CSV.read(dir * "nist_ioniz.txt", header=header, delim=" ",
                      ignorerepeated=true, silencewarnings=true)
        col_to_float!(df, names(df)[3:end]...)
    end
    return df
end

const dfi = tabulate_ionization(nist=true)

"""

"""
function tabulate_VALIIIc(dir::String=datdir)
    @assert isdir(dir)

    # read it in
    header = ["h", "m", "τ_500", "T", "V", "n_H", "n_e", "Ptot", "Pg_Ptot", "ρ"]
    df = CSV.read(dir*"VALIIIC.txt", header=header, comment="#",
                  delim=" ", ignorerepeated=true)
    col_to_float!(df, names(df)...)
    return df
end

const dfv = tabulate_VALIIIc()

function interp_valIIIc(;npoints::Int=1000)
    # make hnew
    ind = 10
    hnew = range(dfv.h[ind], dfv.h[end], length=npoints)
    return interp_valIIIc(hnew)
end


function interp_valIIIc(hnew::AA{T,1}) where T<:Real
    # pre-make df
    dfv_new = DataFrame()
    dfv_new.h = hnew

    # interpolate on each variable in old df
    for key in names(dfv)
        # pass on h
        key == :h && continue

        # interpolate
        spl = Spline1D(reverse(dfv[!, :h]), reverse(dfv[!, key]), k=1, bc="error")
        dfv_new[!, key] = spl(dfv_new.h)
    end
    return dfv_new
end

# function interp_valIIIc(;npoints::Int=1000)
#     # make hnew
#     ind = 2
#     τnew = logspace(dfv.τ_500[ind], dfv.τ_500[end]-0.001, length=npoints)
#     return interp_valIIIc(τnew)
# end


# function interp_valIIIc(τnew::AA{T,1}) where T<:Real
#     # pre-make df
#     dfv_new = DataFrame()
#     dfv_new.τ_500 = τnew

#     # interpolate on each variable in old df
#     for key in names(dfv)
#         # pass on h
#         key == :τ_500 && continue

#         # interpolate
#         spl = Spline1D(dfv.τ_500, dfv[!, key], k=1, bc="error")
#         dfv_new[!, key] = spl(dfv_new.τ_500)
#     end
#     return dfv_new
# end
