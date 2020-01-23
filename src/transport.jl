Sν_lin(τ::T, an...) where T<:Real = an[1] + an[2] * τ
Sν_quad(τ::T, an...) where T<:Real = an[1] + an[2] * τ + an[3] * τ^2
Cν(τ::T, Sν::Function, an...) where T<:Real = Sν(τ, an...) * exp(-τ)
