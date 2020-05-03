"""

"""
function tabulate_departure(dir::String=datdir)
    @assert isdir(dir)

    header = ["h", "b"]
    df = CSV.read(dir * "3sdep_extrap.txt", header=header, delim=" ", ignorerepeated=true)
    col_to_float!(df, names(df)...)
    return df
end

const dfd = tabulate_departure()

"""

"""
function tabulate_source(dir::String=datdir)
    @assert isdir(dir)

    header = ["h", "S"]
    df = CSV.read(dir * "3sdep_extrap.txt", header=header, delim=" ", ignorerepeated=true)
    col_to_float!(df, names(df)...)
    return df
end

const dfs = tabulate_source()

"""

"""
function interp_nlte(df::DataFrame; npoints::Int=1000)
    # make hnew
    ind = 1
    hnew = range(df.h[ind], df.h[end], length=npoints)
    return interp_nlte(df, hnew)
end


function interp_nlte(df::DataFrame, hnew::AA{T,1}) where T<:Real
    # pre-make df
    df_new = DataFrame()
    df_new.h = hnew

    # interpolate on each variable in old df
    for key in names(df)
        # pass on h
        key == :h && continue

        # interpolate
        spl = Spline1D(reverse(df[!, :h]), reverse(df[!, key]), k=1, bc="error")
        df_new[!, key] = spl(df_new.h)
    end
    return df_new
end
