import Bridge.expint

"""
Modified from the source at:
https://github.com/mschauer/Bridge.jl/blob/master/src/expint.jl

Fixed issue where expint(n, 0.0) for n>1 returns NaN
"""
function expint(n::Integer, z)
    if n == 1
        return expint(z)
    elseif n < 1
        # backwards recurrence from E₀ = e⁻ᶻ/z
        zinv = inv(z)
        e⁻ᶻ = exp(-z)
        Eᵢ = zinv * e⁻ᶻ
        for i in 1:-n
            Eᵢ = zinv * (e⁻ᶻ + i * Eᵢ)
        end
        return Eᵢ
    elseif n > 1
        # forwards recurrence from E₁
        e⁻ᶻ = exp(-z)
        Eᵢ = expint(z)
        Eᵢ *= !isinf(Eᵢ)
        for i in 2:n
            Eᵢ = (e⁻ᶻ - z*Eᵢ) / (i - 1)
        end
        return Eᵢ
    end
end
