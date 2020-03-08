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
Ts = [4310.0, 4325.0, 4345.0, 4370.0, 4405.0, 4445.0, 4488.0, 4524.0, 4561.0, 4608.0, 4660.0, 4720.0]
Pg = [2.87, 3.03, 3.17, 3.29, 3.41, 3.52, 3.64, 3.75, 3.86, 3.97, 4.08, 4.19]

# now calculate
Pe = log10.(SA.calc_Pe.(Ts, exp10.(Pg)))

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.plot(Ts, Pe)
ax2.set_xlim(ax1.get_xlim())

# find matching values
spl = Spline1D(Ts, Pg)
Pgs = spl(ax1.get_xticks())
ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
ax2.set_xticklabels(Pgs)
# ax2.set_xticklabels(string.(Pg))
