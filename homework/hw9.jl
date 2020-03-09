using Pkg; Pkg.activate(".")
using Revise
using Dierckx
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# get the electron pressure
Pg = SA.dfv.Pg_Ptot .* SA.dfv.Ptot
Pe0 = SA.dfv.Ptot .- Pg
Pe = SA.calc_Pe.(SA.dfv.T, SA.dfv.Pg_Ptot .* SA.dfv.Ptot, Pe0)

# evaluate Saha for hydrogen & get proton number density
n_p = SA.ΦT.(SA.dfv.T, "H") .* SA.dfv.n_H ./ Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(SA.dfv.h, log10.(SA.dfv.n_e), "k-", label=L"n_{\rm e}")
ax1.plot(SA.dfv.h, log10.(n_p), "k--", label=L"n_{\rm p}")
ax1.set_xlim(0,800)
ax1.set_ylim(8,14)
ax1.invert_xaxis()
ax1.set_xlabel(L"{\rm Height\ (km)}")
ax1.set_ylabel(L"\log(n)")
ax1.legend()
fig.savefig(outdir*"hw9_np_ne.pdf")
plt.clf(); plt.close()
