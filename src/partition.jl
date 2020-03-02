"""

"""
function temp_to_theta(temp::T) where T<:Real
    return 5040.0/(temp)
end

function theta_to_temp(θ::T) where T<:Real
    return 5040.0/θ
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
