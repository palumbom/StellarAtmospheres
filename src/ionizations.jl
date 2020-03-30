"""

"""
function I(species::String, ion::Symbol)
    return χ(species::String, ion::Symbol)
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

Gray Eq. 8.3
"""
function χ(n::Int, species::String)
    return χ(species, :First) - ((h*R*c)/n^2) * eV
end
