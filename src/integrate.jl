"""
    trap_int(x,y)

Compute the trapezoidal integral for pre-evaluated function y=f(x).
"""
function trap_int(x::AbstractArray{T,1}, y::AbstractArray{T,1}) where T<:Real
    dx = x[2:end] .- x[1:end-1]
    num = y[2:end] .+ y[1:end-1]
    return sum(dx .* num ./ 2.0)
end

"""
    trap_int(f, ab; ntrap=NaN, err=false)

Compute the trapezoidal integral for function f(x) on the range (a,b) using
ntrap steps. If err is true, calculate the error via asymptotic estimator
(see wikipedia for summary).
"""
function trap_int(f::Function, ab::Tuple{T,T}; ntrap::Int=NaN, logx::Bool=false, err::Bool=false) where T<:Real
    @assert ntrap > 1
    @assert ab[2] > ab[1]

    # get x-array
    if logx
        x = range(log(ab[1]), log(ab[2]), length=ntrap)
        y = exp.(x) .* f.(exp.(x))
    else
        x = range(ab[1], ab[2], length=ntrap)
        y = f.(x)
    end

    # evaluate integral
    int = sum(diff(x) .* asum(y))/2.0

    # evaluate error
    if err
        db = (y[end] - y[end-1])/(x[end] - x[end-1])
        da = (y[2] - y[1])/(x[2] - x[1])
        err = -((ab[2] - ab[1])^3 / (12*ntrap^2)) * (db - da)
        return int, err
    end
    return int
end
