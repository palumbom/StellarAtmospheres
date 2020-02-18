import QuadGK.quadgk

"""
    trap_int(x,y)

Compute the trapezoidal integral for pre-evaluated function y=f(x).
"""
function trap_int(x::AA{T,1}, y::AA{T,1}) where T<:Real
    return sum(diff(x).*asum(y))/2.0
end

"""
    trap_int(f, ab; ntrap=NaN, logx=false, err=false)

Compute the trapezoidal integral for function f(x) on the range (a,b) using
ntrap steps. If logx is true, integrate on logarithmically-spaced points as
int f(x)xd(lnx). If err is true, calculate the error via asymptotic estimator
(see wikipedia for summary).
"""
function trap_int(f::Function, ab::Tuple{T,T}; ntrap::Int=1, logx::Bool=false, err::Bool=false) where T<:Real
    @assert ntrap > 1
    @assert ab[2] > ab[1]

    # get x-array
    x = range(ab[1], ab[2], length=ntrap)
    if logx
        y = x .* f.(x)
        x = range(log(ab[1]), log(ab[2]), length=ntrap)
    else
        y = f.(x)
    end

    # evaluate integral
    if !err
        return sum(diff(x).*asum(y))/2.0
    else
        db = (y[end] - y[end-1])/(x[end] - x[end-1])
        da = (y[2] - y[1])/(x[2] - x[1])
        err = -((ab[2] - ab[1])^3 / (12*ntrap^2)) * (db - da)
        return sum(diff(x).*asum(y))/2.0, err

    end
end

