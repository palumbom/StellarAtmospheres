function trap_int(x::AbstractArray{T,1}, y::AbstractArray{T,1}) where T<:Real
    dx = x[2:end] .- x[1:end-1]
    num = y[2:end] .+ y[1:end-1]
    return sum(dx .* num ./ 2.0)
end
