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

# some temps and gas pressures
Ts = [4310.0, 4325.0, 4345.0, 4370.0, 4405.0, 4445.0, 4488.0, 4524.0, 4561.0,
      4608.0, 4660.0, 4720.0, 4800.0, 4878.0, 4995.0, 5132.0, 5294.0, 5490.0,
      5733.0, 6043.0, 6429.0, 6904.0, 7467.0, 7962.0, 8358.0, 8630.0, 8811.0]
Pg = [2.87, 3.03, 3.17, 3.29, 3.41, 3.52, 3.64, 3.75, 3.86, 3.97, 4.08, 4.19,
      4.30, 4.41, 4.52, 4.63, 4.74, 4.85, 4.95, 5.03, 5.10, 5.15, 5.18, 5.21,
      5.23, 5.26, 5.29]

# now calculate
Pe = log10.(SA.calc_Pe.(Ts, exp10.(Pg)))

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(Ts, Pe)
spl = Spline1D(Ts, Pg)
Pgs = spl(ax1.get_xticks())
ax2.set_xlim(ax1.get_xlim())
ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
ax2.set_xticklabels(round.(Pgs, digits=2))
ax2.set_xlabel(L"\log P_{\rm g}\ ({\rm dyne\ cm}^{-2})")
ax1.set_xlabel(L"T\ {\rm (K)}")
ax1.set_ylabel(L"\log P_{\rm e}\ ({\rm dyne\ cm}^{-2})")
fig.savefig(outdir*"hw8_Pe.pdf")
plt.clf(); plt.close()
