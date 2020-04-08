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
θ = range(0.2, 2.0, step=0.2)
T = SA.theta_to_temp.(θ)

# interpolation grid
Ts = range(2520.0, 25200.0, length=100)

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(T, SA.row_for_species(df, "Si"), c="k", label="Silicon")
ax1.scatter(T, SA.row_for_species(df, "Fe"), c="b", label="Iron")
ax1.scatter(T, SA.row_for_species(df, "Mg"), c="r", label="Magnesium")
ax1.plot(Ts, log10.(SA.calc_partition.(Ts, "Si")), "k--")
ax1.plot(Ts, log10.(SA.calc_partition.(Ts, "Fe")), "b--")
ax1.plot(Ts, log10.(SA.calc_partition.(Ts, "Mg")), "r--")
ax1.set_xlabel(L"{\rm Temperature\ (K)}")
ax1.set_ylabel(L"\log u(T)")
ax1.set_xlim(1500, 26500)
ax1.legend()
fig.savefig(outdir*"hw6_partit.pdf")
plt.clf(); plt.close()
