function deriv(f::Function, x::AA{T,1}) where T<:Real
    dx = x[2] - x[1]
    return (f.(x .+ dx) .- f.(x .- dx))./(2*dx)
end

function deriv2(f::Function, x::AA{T,1}) where T<:Real
    dx = x[2] - x[1]
    num = f.(x .+ dx) .- 2 .* f.(x) .+ f.(x .- dx)
    return num./(dx^2)
end
