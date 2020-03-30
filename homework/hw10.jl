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

# get VALIIIC atmosphere
dfv = SA.dfv

# get stuff at τ=1 surface
ind = 49
temp = dfv.T[ind]
Pg = dfv.Pg_Ptot[ind] .* dfv.Ptot[ind]
Pe1 = SA.P_from_nkt.(dfv.n_e[ind], temp)
Pe2 = SA.calc_Pe.(temp, Pg, Pe1)
ξ = dfv.V[49]

# make LineParams object instances for NaD lines
m = 3.817541e-23
NaD1 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1], m=m, gu=4, gl=2, logC4=-15.33)
NaD2 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1], m=m, gu=2, gl=2, logC4=-15.17)

# get Voigt profile for line 1
waves = range(5880.0, 5900.0, length=1000)
vprof = SA.calc_α(5890.0, Pe2, Pg, temp, ξ, NaD1)
