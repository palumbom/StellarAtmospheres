function tabulate_abundances(dir::String=datdir)
    @assert isdir(dir)

    # read it in
    header = ["Z", "Element", "Weight", "A", "logA", "logA12"]
    df = CSV.read(dir*"SolarAbundance.txt", header=header, skipto=2)
    col_to_float!(df, names(df)[4:end]...)
    return df
end

const dfa = tabulate_abundances()
