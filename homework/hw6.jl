using Pkg; Pkg.activate(".")
using Revise
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# read in the partition and plot example
df = SA.tabulate_partition()

# get species
s = "Si"
θ = range(0.2, 2.0, step=0.2)
T = SA.theta_to_temp.(θ)
P = SA.row_for_species(df, s)

# interpolation grid
Ts = range(2520.0, 25200.0, length=100)

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(T, P, c="k", label="Data")
ax1.plot(Ts, SA.calc_partition(Ts, df, s), "k--", label="Interpolation")
ax1.set_xlabel(L"{\rm Temperature\ (K)}")
ax1.set_ylabel(L"{\rm Partition\ Function}")
ax1.legend()
fig.savefig(outdir*"hw6_partit.pdf")
