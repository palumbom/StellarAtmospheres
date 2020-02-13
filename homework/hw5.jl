using Pkg; Pkg.activate(".")
using Revise
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# parameters for flux calc
ν̃ = range(0.01, 12.0, length=5000)
ν = ν̃2ν.(ν̃.*1e4)
T = 8700.0
τs = (1e-10, 20.0)

# calculate the emergent flux
F8777 = ℱν₀(SνPlanck, τs, Teff=T, ν=ν, ntrap=500)

# plot the emergent flux
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(ν̃, F8777, label="Integral")
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"{\rm Emergent Flux}\ {\rm erg/s/cm}^2{\rm Hz/ster}")
ax1.set_ylim(1e-7, 1e-3)
ax1.set_yscale("log")
ax1.legend()
ax1.annotate(L"T_{\rm eff} = 8777\ {\rm K}", (0.1, 5e-4), xycoords="data")
fig.savefig(outdir * "hw5_f8777.pdf")
plt.clf(); plt.close()
