# define composite type to hold line characteristics
struct LineParams{T<:AF}
    element::String # name of element
    n::Integer      # energy leve of lower level
    λ₀::T           # rest frame line center
    ν₀::T           # line center but ~frequency~
    A::AA{T,1}      # Einstein A coefficients
    m::T            # mass of atomic in cgs
    gu::Integer     # statistical weight of upper level
    gl::Integer     # statistical weight of lower level
    logC4::T        # stark interaction constant
    logC6::T        # van der waals interaction constant
end

# make an outer constructor so you can't mess up arg order
function LineParams(;element="", n=NaN, λ₀=NaN, A=[NaN], m=NaN, gu=NaN, gl=NaN, logC4=NaN)
    ν₀ = λ2ν(λ₀)
    logC6 = calc_logC6(element, n, λ₀)
    return LineParams(element, n, λ₀, ν₀, A, m, gu, gl, logC4, logC6)
end

"""
    calc_logC6()

Calculate the log of C6 for the Van der Waals dispersion. Gray Eq. 11.30.
"""
function calc_logC6(element::String, n::Int, λ₀::T) where T<:AF
    # get ionization/excitation stuff
    I1 = I(element, :First)
    χ1 = χ(n, "Na")
    χ1 = 0.0        # gives -31.7 which agrees with Gray for unknown reasons
    χλ = 1.2398e4/λ₀

    # calculate terms and return
    term1 = (I1 - χ1 - χλ)^2
    term2 = (I1 - χ1)^2
    return log10(0.3e-30 * (one(T)/term1 - one(T)/term2))
end
