using Pkg; Pkg.activate(".")
using Dierckx
using BenchmarkTools
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# interpolate departure coefficients and source functions
npoints = 1000
depc_df = SA.interp_nlte(SA.dfd, npoints=npoints)
source_df = SA.interp_nlte(SA.dfs, npoints=npoints)

# actuall pull out the height, b, S arrays
h = depc_df.h
b = depc_df.b
S = source_df.S
