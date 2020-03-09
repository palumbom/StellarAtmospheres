function estimate_Pe(temp::T, Pg::T) where T<:AF
    t1 = 0.5 * Pg
    t2 = ΦT(temp, "H")
    return minimum([t1, t2])
end

function calc_Pe(temp::T, Pg::T; atol=1e-4) where T<:AF
    # get initial pressure estimate
    Pe0 = estimate_Pe(temp, Pg)

    # pre-calc some stuff
    Φj = ΦT.(temp, dfa.Element)
    Aj = dfa.A

    # set initial Pe and iterate
    Pe = Pe0
    i = 0
    while true
        i += 1
        num = sum(filter(!isnan, Aj.*(Φj./Pe)./(one(T).+(Φj./Pe))))
        den = sum(filter(!isnan, Aj.*(one(T).+(Φj./Pe)./(one(T) .+ (Φj./Pe)))))

        # return if converged
        if abs(Pe - Pg * num/den) <= atol
            # @show i
            return Pg * num/den
        end

        # else update and interate
        Pe = Pg * num/den
    end
end

function calc_Pe(temp::T, Pg::T, Pe0::T; atol=1e-4) where T<:AF
    @assert !isnan(Pe0)

    # pre-calc some stuff
    Φj = ΦT.(temp, dfa.Element)
    Aj = dfa.A

    # set initial Pe and iterate
    Pe = Pe0
    i = 0
    while true
        i += 1
        num = sum(filter(!isnan, Aj.*(Φj./Pe)./(one(T).+(Φj./Pe))))
        den = sum(filter(!isnan, Aj.*(one(T).+(Φj./Pe)./(one(T) .+ (Φj./Pe)))))

        # return if converged
        if abs(Pe - Pg * num/den) <= atol
            # @show i
            return Pg * num/den
        end

        # else update and interate
        Pe = Pg * num/den
    end
end

# function
