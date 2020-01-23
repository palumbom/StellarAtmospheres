using Pkg; Pkg.activate(".")
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib;
AA = AbstractArray

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# define some source functions
Sν_lin(τ::T, an...) where T<:Real = an[1] .+ an[2] .* τ
Sν_quad(τ::T, an...) where T<:Real = an[1] .+ an[2] .* τ .+ an[3] .* τ.^2

# evaulate the source functions
τ  = range(0.0, 100.0, length=1000)
an = [1.0, 1.0, 1.0]
S1 = Sν_lin.(τ, an...)
S2 = Sν_quad.(τ, an...)

# plot the source functions
