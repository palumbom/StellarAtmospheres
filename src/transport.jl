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
coefficients a_n.
Coefficients should be passed as multiple arguments or a splatted array.
"""
Cν(τ::T, an::T...) where T<:Real = Sν(τ, an...) * exp(-τ)
Cν(τ::T, μ::T, an::T...) where T<:Real = Sν(τ, an...) * exp(-τ/μ) / μ

"""
    Iν₀(μ, an...)

Compute the emergent intensity at μ by integration.
Coefficients should be passed as multiple arguments or a splatted array.
"""
function Iν₀(μ::T, τs::Tuple{T,T}, an::T...; ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    return trap_int(x -> Cν(x, μ, an...), τs, ntrap=ntrap, logx=true)
end

"""
    IνEB(μ, an...)

Compute the emergent intensity at μ using the Eddington-Barbier approximation.
Coefficients should be passed as multiple arguments or a splatted array.
"""
function IνEB(μ::T, an::T...) where T<:Real
    return Sν(μ, an...)
end

"""
    ℱν₀(τs::Tuple, an...; ntrap=NaN)

Compute the surface flux (as defined in Rutten) by integrating over τs given
source function Sν with coefficients a_n. Use a trapezoidal integrator with
ntrap trapezoids and logarithmically-spaced gridpoints. Source function
coefficients should be passed as multiple arguments or a splatted array.
"""
function ℱν₀(τs::Tuple{T,T}, an::T...; ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    f = x -> Sν(x, an...) * expint(2, x)
    return (2.0 * π) * trap_int(f, τs, ntrap=ntrap, logx=true)
end

"""
    Hν₀(τs::Tuple, an...; ntrap=NaN)

Compute the first moment of intensity (as defined in Rutten) by
integrating over τs given source function Sν with coefficients a_n.
Use a trapezoidal integrator with ntrap trapezoids and logarithmically-spaced
gridpoints. Source function coefficients should be passed as multiple
arguments or a splatted array.
"""
function Hν₀(τs::Tuple{T,T}, an::T...; ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    return ℱν₀(τs, an..., ntrap=ntrap)/(4.0 * π)
end

"""
    Hν₀(an...; ntrap=NaN, EB=true)

Compute the first moment of intensity (as defined in Rutten) by directly
integrating the emergent intensity over μ given source function Sν
with coefficients a_n. Use a trapezoidal integrator with ntrap trapezoids.
Source function coefficients should be passed as multiple arguments or a
splatted array.
"""
function Hν₀(an::T...; ntrap::Int=NaN, EB::Bool=true) where T<:Real
    @assert !isnan(ntrap)
    μs = (1e-10, 1.0)
    if EB
        return 0.5 * trap_int(x -> x*IνEB(x, an...), μs, ntrap=ntrap)
    else
        τs = (1e-10, 1000.0)
        return 0.5 * trap_int(x -> x*Iν₀(x, τs, an..., ntrap=ntrap), μs, ntrap=ntrap)
    end
end

"""
    Hν₀(an...)

Compute the first moment of intensity (as defined in Rutten) via the
Eddington-Barbier approximation for Hν(0). Source function coefficients
should be passed as multiple arguments or a splatted array.
"""
function HνEB(an::T...) where T<:Real
    return 0.25 * Sν((2.0/3.0), an...)
end
