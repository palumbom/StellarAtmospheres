using Pkg; Pkg.activate(".")
using Revise
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# calculate opacities
λs = range(1000.0, 20000.0, length=1000)

# do Gray Fig. 8.5b first
T = 6429.0
Pe = exp10(1.77)
κHbf = SA.κ_H_bf.(λs, T, Pe)
κHff = SA.κ_H_ff.(λs, T, Pe)
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe)
κHmff = SA.κ_H_minus_ff.(λs, T, Pe)
κT = SA.κ_tot.(λs, T, Pe)

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label="Total")
ax1.plot(λs, κHbf, "b-.", label="H Bound-Free")
ax1.plot(λs, κHff, "b:", label="H Free-Free")
ax1.plot(λs, κHmbf, "r-.", label="H- Bound-Free")
ax1.plot(λs, κHmff, "r:", label="H- Free-Free")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e\ {\rm Units\ tbd}")
ax1.legend()
fig.savefig(outdir * "hw7_85b.pdf")

# do Gray Fig. 8.5c now
T = 6429.0
Pe = exp10(1.77)
κT = SA.κ_tot.(λs, T, Pe)
