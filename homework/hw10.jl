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

# wavelengths
NaD1 = 5890.0
NaD2 = 5896.0

# make LineParams object instances for NaD lines
m = 3.817541e-23
# NaD1 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1], m=m, gu=4, gl=2, logC4=-15.33)

"""
NaD[0]['wavelength'] = 5890*u.AA
NaD[1]['wavelength'] = 5896*u.AA

NaD[0]['A21']=1E8*6.16E-1/u.s
NaD[1]['A21']=1E8*6.14E-1/u.s

NaD[0]['gu'] = 4
NaD[1]['gu'] = 2

NaD[0]['gl'] = 2
NaD[1]['gl'] = 2
"""
