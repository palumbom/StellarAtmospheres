"""

"""
function temp_to_theta(temp::T) where T<:Real
    return 5039.78/(temp)
end

function theta_to_temp(θ::T) where T<:Real
    return 5039.78/θ
end

"""
Calculate partition functions from Gray using his Table D.2 to interpolate
"""
function calc_partition(temp::T, species::String) where T<:AF
    @assert temp >= zero(T)

    # check for simple species
    if species == "H-"
        return zero(T)
    elseif species == "HII"
        return zero(T)
    end

    # now check for species in table
    if !(species in dfp.Species)
        return zero(T)
    end

    # theta data
    θs = range(0.2, 2.0, step=0.2)
    Ps = row_for_species(dfp, species)

    # now do the interpolation
    ninds = .!isnan.(Ps)
    spl = Spline1D(θs[ninds], Ps[ninds], k=1, bc="nearest") # linear interp log debug
    return spl(temp_to_theta(temp))
end
