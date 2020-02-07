import Bridge.expint

function expint_log(n::Int, u::T) where T<:Real
    return log(10) * exp(u) * expint(n, exp(u))
end
