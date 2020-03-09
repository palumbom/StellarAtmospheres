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

# shortcut to data
dfv = SA.dfv
dfa = SA.dfa

# get the electron pressure from ideal gas law
Pe1 = SA.P_from_nkt.(dfv.n_e, dfv.T)

# also try electron pressure from iterative solving
Pg = dfv.Pg_Ptot .* dfv.Ptot
Pe0 = dfv.Ptot .- Pg
Pe2 = SA.calc_Pe.(dfv.T, dfv.Pg_Ptot .* dfv.Ptot, Pe0)

# evaluate Saha for hydrogen & get proton number density
n_p1 = SA.ΦT.(dfv.T, "H") .* dfv.n_H ./ Pe1
n_p2 = SA.ΦT.(dfv.T, "H") .* dfv.n_H ./ Pe2

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(dfv.h, log10.(dfv.n_e), "k-", label=L"n_{\rm e}")
ax1.plot(dfv.h, log10.(n_p1), "k--", label=L"n_{\rm p}")
# ax1.plot(dfv.h, log10.(n_p2), "k-.", label=L"n_{\rm p}")
ax1.set_xlim(0,800)
ax1.set_ylim(8,14)
ax1.invert_xaxis()
ax1.set_xlabel(L"{\rm Height\ (km)}")
ax1.set_ylabel(L"\log(n)")
ax1.legend()
fig.savefig(outdir*"hw9_np_ne.pdf")
plt.clf(); plt.close()

# get number density of electrons from each element
n_e_H = SA.ΦT.(dfv.T, "H").* dfv.n_H .* SA.abundance_for_element("H") ./ Pe2
n_e_Fe = SA.ΦT.(dfv.T, "Fe").* dfv.n_H .* SA.abundance_for_element("Fe") ./ Pe2
n_e_Mg = SA.ΦT.(dfv.T, "Mg").* dfv.n_H .* SA.abundance_for_element("Mg") ./ Pe2
n_e_Si = SA.ΦT.(dfv.T, "Si").* dfv.n_H .* SA.abundance_for_element("Si") ./ Pe2

# total electron density is sum of each
n_e_tot = n_e_H .+ n_e_Fe .+ n_e_Mg .+ n_e_Si

# fraction
f_H = n_e_H./n_e_tot
f_Fe = n_e_Fe./n_e_tot
f_Mg = n_e_Mg./n_e_tot
f_Si = n_e_Si./n_e_tot

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(dfv.h, f_H, "k-", label="H")
ax1.plot(dfv.h, f_Fe, "k-.", label="Fe")
ax1.plot(dfv.h, f_Mg, "k--", label="Mg")
ax1.plot(dfv.h, f_Si, "k:", label="Si")
ax1.set_xlim(0,800)
ax1.invert_xaxis()
ax1.set_xlabel(L"{\rm Height\ (km)}")
ax1.set_ylabel(L"{\rm Contribution\ to\ } n_{\rm e}")
ax1.legend()
fig.savefig(outdir*"hw9_f_ne.pdf")
plt.clf(); plt.close()

# plot the electron pressures
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(dfv.h, Pe1, "k--", label="Ideal Gas")
ax1.plot(dfv.h, Pe2, "k-", label="Iterative")
ax1.set_xlim(0,800)
ax1.invert_xaxis()
ax1.set_yscale("log")
ax1.set_xlabel(L"{\rm Height\ (km)}")
ax1.set_ylabel(L"P_e\ ({\rm dyne\ cm}^{-2})")
ax1.legend()
fig.savefig(outdir*"hw9_Pe.pdf")
plt.clf(); plt.close()

# parameters for total opacity
T = 6420.0
Pe = 57.0
Pg = 1.13e5
λ = 5000.0

# stimulated emission and the conversion term
SE = (1.0 .- exp10.(-1.2398e4./λ.*SA.temp_to_theta(T)))
PT = (1.0 + SA.ΦT(T, "H")/Pe)
pg = (mH * sum(filter(!isnan, dfa.A .* dfa.Weight)))

# neutral fraction
frac = SA.neutral_fraction(T, Pe, "H")

# sources of opacity
κcont = SA.κ_tot(λ, T, Pe, Pg)
κHmbf = SA.κ_H_minus_bf(λ, T, Pe) .* SE ./ PT ./ pg
κHmff = SA.κ_H_minus_ff(λ, T, Pe) ./ PT ./ pg
κHbf = SA.κ_H_bf(λ, T, Pe) .* SE ./ PT ./ pg
κHff = SA.κ_H_ff(λ, T, Pe) .* SE ./ PT ./ pg
κe = SA.κ_e(Pe, Pg) ./ pg

# calculate continuum opacity as function of Pe
κconts = SA.κ_tot.(λ, dfv.T, Pe2, dfv.Pg_Ptot .* dfv.Ptot)

# now plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(dfv.τ_500, κconts, "k-")
ax1.set_yscale("log")
ax1.set_xlabel(L"\tau_{500}")
ax1.set_ylabel(L"\kappa_{500}\ ({\rm cm}^2\ {\rm g}^{-1})")
fig.savefig(outdir*"hw9_kappa.pdf")
plt.clf(); plt.close()

# verify opacity
g = exp10(4.4377)
dτ = diff(dfv.τ_500)
dP = diff(dfv.Ptot)
