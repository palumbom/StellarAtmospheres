"""
    Sν(τ, an...)

Compute the source function at τ with coefficients a_n.
Coefficients should be passed as multiple arguments or a splatted array.
"""
function Sν(τ::T, an::Vararg{T,N}) where {T<:Real, N}
    return sum([an[n] * τ^(n-1) for n in 1:N])
end

"""
    Cν(τ, an...)

Compute the contribution function at τ with the specified source function and
coefficients an.
Coefficients should be passed as multiple arguments or a splatted array.
"""
Cν(τ::T, Sν::Function, an::T...) where T<:Real = Sν(τ, an...) * exp(-τ)

"""
    ℱν₀(τs::Tuple, Sν, an...; ntrap=NaN)

Compute the surface flux (as defined in Rutten) by integrating over τs given
source function Sν with coefficients an. Use a trapezoidal integrator with
ntrap trapezoids and logarithmically-spaced gridpoints. Source function
coefficients should be passed as multiple arguments or a splatted array.
"""
function ℱν₀(τs::Tuple{T,T}, Sν::Function, an::T...; ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    f = x -> Sν(x, an...) * expint(2, x)
    return 2π * trap_int(f, (τs[1], τs[end]), ntrap=ntrap, logx=true)
end

function Hν₀(τs::Tuple{T,T}, Sν::Function, an::T...; ntrap::Int=NaN) where T<:Real
    return ℱν₀(τs, Sν, an..., ntrap=ntrap)/(4π)
end
