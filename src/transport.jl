"""
    SνPoly(τ; an=[NaN])

Compute the source function at τ as a polynomial expansion with
coefficients a_n.
"""
function SνPoly(τ::T; an::AA{T,1}=[NaN], kwargs...) where T<:Real
    @assert sum(isnan.(an)) == 0
    return sum([an[n] * τ^(n-1) for n in 1:length(an)])
end

"""
    SνPlanck(τ, ν, Teff=NaN)

Compute the source function at τ and ν as the Planck Function with Teff.
"""

function SνPlanck(τ::T; Teff::T=NaN, ν::AA{T,1}=[NaN], kwargs...) where T<:Real
    @assert !isnan(Teff)
    @assert sum(isnan.(ν)) == 0
    Ts = Tτ(τ, Teff=Teff)
    return Bν.(ν, Ts)
end

"""
    Cν(τ, an...; μ=1.0, Teff, ν)

Compute the contribution function at optical depth τ and and disk position
μ with the specified source function. For a polynomial Sν, coefficients
an... must be passed. For Planck Function, an effective temperature and an
array of frequencies must be passed.
"""
function Cν(Sν::Function, τ::T; μ::T=1.0, Teff::T=NaN, an::AA{T,1}=[NaN], ν::AA{T,1}=[NaN]) where T<:Real
    return Sν(τ, an=an, ν=ν, Teff=Teff) * exp(-τ/μ)/μ
end

"""
    Iν₀(μ, τs, an...)

Compute the emergent intensity at μ by integration.
Coefficients should be passed as multiple arguments or a splatted array.
"""
function Iν₀(Sν::Function, μ::T, τs::Tuple{T,T}; Teff::T=NaN, an::AA{T,1}=[NaN], ν::AA{T,1}=[NaN], ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    f = x -> Cν(Sν, x, an=an, Teff=Teff, ν=ν, μ=μ)
    return trap_int(f, τs, ntrap=ntrap, logx=true)
end

"""
    IνEB(μ, an...)

Compute the emergent intensity at μ using the Eddington-Barbier approximation.
Coefficients should be passed as multiple arguments or a splatted array.
"""
function IνEB(Sν::Function, μ::T; Teff::T=NaN, ν::AA{T,1}=[NaN], an::AA{T,1}=[NaN]) where T<:Real
    return Sν(μ, an=an, ν=ν, Teff=Teff)
end

"""
    ℱν₀(τs::Tuple, an...; ntrap=NaN)

Compute the surface flux (as defined in Rutten) by integrating over τs given
source function Sν with coefficients a_n. Use a trapezoidal integrator with
ntrap trapezoids and logarithmically-spaced gridpoints. Source function
coefficients should be passed as multiple arguments or a splatted array.
"""
function ℱν₀(Sν::Function, τs::Tuple{T,T}; Teff::T=NaN, an::AA{T,1}=[NaN], ν::AA{T,1}=[NaN], ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    f1 = x -> Sν(x, an=an, ν=ν, Teff=Teff) .* expint(2, x)
    return (2.0 * π) .* trap_int(f1, τs, ntrap=ntrap, logx=false)
end


function ℱντ(Sν::Function, τ::T, τs::Tuple{T,T}; Teff::T=NaN, an::AA{T,1}=[NaN], ν::AA{T,1}=[NaN], ntrap::Int=NaN) where T<:Real
    @assert τs[1] <= τ <= τs[2]
    @assert !isnan(ntrap)

    f1 = x-> Sν(x, an=an, ν=ν, Teff=Teff) .* expint(2, x - τ)
    f2 = x-> Sν(x, an=an, ν=ν, Teff=Teff) .* expint(2, τ - x)
    ugtν = trap_int(f1, (τ, τs[2]), ntrap=ntrap, logx=true)
    dgtν = trap_int(f2, (τs[1], τ), ntrap=ntrap, logx=false)
    return (2.0 * π) .* (ugtν .- dgtν)
end

function ℱτ(Sν::Function, τ::T, τs::Tuple{T,T}; Teff::T=NaN, an::AA{T,1}=[NaN], ν::AA{T,1}=[NaN], ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    return trap_int(ν, ℱντ(Sν, τ, τs, an=an, Teff=Teff, ν=ν, ntrap=ntrap))
end

"""
    Hν₀(τs::Tuple, an...; ntrap=NaN)

Compute the first moment of intensity (as defined in Rutten) by
integrating over τs given source function Sν with coefficients a_n.
Use a trapezoidal integrator with ntrap trapezoids and logarithmically-spaced
gridpoints. Source function coefficients should be passed as multiple
arguments or a splatted array.
"""
function Hν₀(Sν::Function, τs::Tuple{T,T}; an::AA{T,1}=[NaN], Teff::T=NaN, ν::AA{T,1}=[NaN], ntrap::Int=NaN) where T<:Real
    @assert !isnan(ntrap)
    return ℱν₀(Sν, τs, an=an, ν=ν, Teff=Teff, ntrap=ntrap)/(4.0 * π)
end

"""
    Hν₀(an...; ntrap=NaN, EB=true)

Compute the first moment of intensity (as defined in Rutten) by directly
integrating the emergent intensity over μ given source function Sν
with coefficients a_n. Use a trapezoidal integrator with ntrap trapezoids.
Source function coefficients should be passed as multiple arguments or a
splatted array.
"""
function Hν₀(Sν::Function; τs::Tuple{T,T}=(1e-10, 1000.0),
             μs::Tuple{T,T}=(1e-10, 1.0), Teff::T=NaN,
             an::AA{T,1}=[NaN], ν::AA{T,1}=[NaN], ntrap::Int=NaN,
             EB::Bool=true, logx::Bool=true) where T<:Real
    @assert !isnan(ntrap)
    @assert μs[1] >= 0.0
    if EB
        return 0.5 * trap_int(x -> x*IνEB(x, an=an), μs, ntrap=ntrap, logx=logx)
    else
        f = x -> x*Iν₀(Sν, x, τs, an=an, Teff=Teff, ν=ν, ntrap=ntrap)
        return 0.5 * trap_int(f, μs, ntrap=ntrap, logx=logx)
    end
end

"""
    Hν₀(an...)

Compute the first moment of intensity (as defined in Rutten) via the
Eddington-Barbier approximation for Hν(0). Source function coefficients
should be passed as multiple arguments or a splatted array.
"""
function HνEB(Sν; Teff::T=NaN, ν::AA{T,1}=[NaN], an::AA{T,1}=[NaN]) where T<:Real
    return 0.25 * Sν((2.0/3.0), an=an, ν=ν, Teff=Teff)
end

function ℱνEB(Sν; Teff::T=NaN, ν::AA{T,1}=[NaN], an::AA{T,1}=[NaN]) where T<:Real
    return (4.0 * π) * HνEB(Sν, an=an, ν=ν, Teff=Teff)
end
