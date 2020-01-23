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
ax1.set_xlim(0,10)
ax1.legend()
fig.savefig(outdir*"hw3_source.pdf")
plt.clf(); plt.close()

# now define the contribution function
Cν(τ::T, Sν::Function, an...) where T<:Real = Sν(τ, an...) * exp(-τ)

# evaluate the contribution function
C1 = Cν.(τ, Sν_lin, an...)
C2 = Cν.(τ, Sν_quad, an...)

# plot the integrand
fig = plt.figure("Contribution Functions")
ax1 = fig.add_subplot(111)
ax1.plot(τ, C1, label="Linear " * L"S_\nu")
ax1.plot(τ, C2, label="Quadratic " * L"S_\nu")
ax1.set_xlabel(L"\tau_\nu")
ax1.set_ylabel(L"S_\nu(\tau_\nu) \exp(\tau_\nu)")
ax1.set_xlim(0,10)
ax1.legend()
fig.savefig(outdir*"hw3_contribution.pdf")
plt.clf(); plt.close()

# integrate out the contribution functions
nt = 10000           # number of trapezoids
τs = (0.0, 100.0)    # integral bounds for τ
I1 = trap_int(x -> Cν(x, Sν_lin, an...), τs, nt)
I2 = trap_int(x -> Cν(x, Sν_quad, an...), τs, nt)

# compare to E-B approx
println(Sν_lin(1.0, an...))
println(Sν_quad(1.0, an...))
