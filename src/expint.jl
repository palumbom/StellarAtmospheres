import Bridge.expint

function expint_log(n::Int, u::T) where T<:Real
    return log(10) * exp10(u) * expint(n, exp10(u))
end
