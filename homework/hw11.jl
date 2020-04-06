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

# get VALIIIc stuff
dfv = SA.dfv
temp = dfv.T
ne = dfv.n_e
Pe = SA.P_from_nkt.(ne, temp)
Pg = dfv.Pg_Ptot .* dfv.Ptot

# wavelength range generator
waves = range(5888.0, 5898.0, length=500)

# get continous opacity
κcont = SA.κ_tot(waves, temp, Pe, Pg)

# get sodium opacity
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")

# make LineParams object instances for NaD lines
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)


# get τ

# get NaD lines
# flux = ℱν₀()
