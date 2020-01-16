# define planck functions
Bν(T::t, ν::t) where t<:Real = (2.0*h*ν^3/c^2) * 1.0/(exp(h*ν/(kB*T)) - 1.0)
Bλ(T::t, λ::t) where t<:Real = (2.0*h*c^2/λ^5) * 1.0/(exp(h*c/(λ*kB*T)) - 1.0)

# methods of planck function given arrays of freq/wave
function Bν(T::t, ν::AbstractArray{t,1}) where t<:Real
    out = zeros(length(ν))
    for i in 1:length(ν)
        out[i] = Bν(T, ν[i])
    end
    return out
end

function Bλ(T::t, λ::AbstractArray{t,1}) where t<:Real
    out = zeros(length(λ))
    for i in 1:length(λ)
        out[i] = Bλ(T, λ[i])
    end
    return out
end
