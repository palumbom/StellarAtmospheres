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
Sν_lin(τ::T, an...) where T<:Real = an[1] + an[2] * τ
Sν_quad(τ::T, an...) where T<:Real = an[1] + an[2] * τ + an[3] * τ^2

# evaulate the source functions
τ  = range(0.0, 10.0, length=1000)
an = [1.0, 1.0, 1.0]
S1 = Sν_lin.(τ, an...)
S2 = Sν_quad.(τ, an...)

# plot the source functions
fig = plt.figure("Source Functions")
ax1 = fig.add_subplot(111)
ax1.plot(τ, S1, label="Linear")
ax1.plot(τ, S2, label="Quadratic")
ax1.set_xlabel(L"\tau_\nu")
ax1.set_ylabel(L"S_\nu(\tau_\nu)")
ax1.legend()
fig.savefig(outdir*"hw3_source.pdf")
plt.clf(); plt.close()

# now define the contribution function integrand
Cν(τ::T, Sν::Function, an...) where T<:Real = Sν(τ, an...) * exp(τ)

# evaluate the contribution functions
C1 = Cν.(τ, Sν_lin, an...)
C2 = Cν.(τ, Sν_quad, an...)
