# constants in cgs units
const c = 2.99792458e10
const h = 6.62606885e-27
const R = 1.0968e5
const e = 4.80320425e-10
const σB = 5.670374419e-5
const kB = 1.38064852e-16
const kBev = 0.0
const eV = 6.2415e11
const α0 = 1.0449e-26
const mp = 1.67262192369e-24
const me = 9.1093837015e-28
const αe = 6.6524587158e-25     # aka thomson cross-section
const mH = 1.6735575e-24
const NA = 6.022e23
const loge = log10(exp(1.0))    # b/c Gray is weird

# wien displacement law (cm and Hz and K)
λmax(temp::T) where T<:Real = 0.290/temp
νmax(temp::T) where T<:Real = 5.88e10 * temp

# convert λ to ν and ν to λ
λ2ν(λ::T) where T<:Real = c/λ
ν2λ(ν::T) where T<:Real = c/ν
λ2ν̃(λ::T) where T<:Real = one(T)/λ
ν̃2λ(ν̃::T) where T<:Real = one(T)/ν̃
ν̃2ν(ν̃::T) where T<:Real = c*ν̃
ν2ν̃(ν::T) where T<:Real = ν/c
