function legendre(n::Int, x::AA{T,1}) where T<:AbstractFloat
    if n == 0
        return ones(eltype(x), length(x))
    elseif n == 1
        return x
    else
        return ((2.0*n-1.0).*x.*legendre(n-1,x) .- (n-1).*legendre(n-2,x))/n
    end
end

function dlegendre(n::Int, x::AA{T,1}) where T<:AbstractFloat
    if n == 0
        return zero(eltype(x))
    elseif n == 1
        return ones(eltype(x), length(x))
    else
        return (n./(x.^2.0.-1.0)) .* (x.*legendre(n,x) - legendre(n-1,x))
    end
end

function roots_legendre(order::Int, tol=1e-20)
    @assert order >= 2

    for i in 1:(order ÷ 2)  # ÷ performs integer division
    end
end
