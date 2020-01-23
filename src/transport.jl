"""
    Sν_lin(τ, an...)

Compute the linear source function at τ with coefficients a0 and a1.
Coefficients should be passed as multiple arguments or a splatted array.
"""
Sν_lin(τ::T, an::T...) where T<:Real = an[1] + an[2] * τ

"""
    Sν_quad(τ, an...)

Compute the quadratic source function at τ with coefficients a0, a1, a2.
Coefficients should be passed as multiple arguments or a splatted array.
"""
Sν_quad(τ::T, an::T...) where T<:Real = an[1] + an[2] * τ + an[3] * τ^2

"""
    Cν(τ, an...)

Compute the contribution function at τ with the specified source function and
coefficients an.
Coefficients should be passed as multiple arguments or a splatted array.
"""
Cν(τ::T, Sν::Function, an::T...) where T<:Real = Sν(τ, an...) * exp(-τ)
