"""

"""
function tabulate_abundances(dir::String=datdir)
    @assert isdir(dir)

    # read it in
    header = ["Z", "Element", "Weight", "A", "logA", "logA12"]
    df = CSV.read(dir*"SolarAbundance.txt", header=header, skipto=2)
    col_to_float!(df, names(df)[4:end]...)
    return df
end

# constant global variable
const dfa = tabulate_abundances()

"""

"""
function abundance_for_element(element::String)
    # find the appropriate row
    ind = findfirst(dfa.Element .== element)
    return dfa[ind, :A]
end

function sum_abundance_weights()
    return sum(filter(!isnan, dfa.A .* dfa.Weight)) / NA
end
