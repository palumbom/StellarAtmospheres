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

# electron pressure
Pg = dfv.Pg_Ptot .* dfv.Ptot
Pe1 = SA.P_from_nkt.(dfv.n_e, dfv.T)
Pe2 = SA.calc_Pe.(dfv.T, dfv.Pg_Ptot .* dfv.Ptot, Pe1)

# evaluate Saha for hydrogen & get proton number density
n_p = SA.species_fraction.(dfv.T, Pe2, "H", ion="First") .* dfv.n_H
n_e_H = n_p

# get the number density from ideal gas law
n_e = Pe2./(SA.kB .* dfv.T)

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(dfv.h, log10.(n_e), "k-", label=L"n_{\rm e}")
ax1.plot(dfv.h, log10.(n_p), "k--", label=L"n_{\rm p}")
ax1.set_xlim(0,800)
ax1.set_ylim(8,14)
ax1.invert_xaxis()
ax1.set_xlabel(L"{\rm Height\ (km)}")
ax1.set_ylabel(L"\log(n)")
ax1.legend()
fig.savefig(outdir*"hw9_np_ne.pdf")
plt.clf(); plt.close()

# do Iron now
n1nt = SA.species_fraction.(dfv.T, Pe2, "Fe", ion="First")
n2nt = SA.species_fraction.(dfv.T, Pe2, "Fe", ion="Second")
n_e_Fe = (n1nt .+ 2.0.*n2nt) .* dfv.n_H .* SA.abundance_for_element("Fe")

# do Silicon now
n1nt = SA.species_fraction.(dfv.T, Pe2, "Si", ion="First")
n2nt = SA.species_fraction.(dfv.T, Pe2, "Si", ion="Second")
n_e_Si = (n1nt .+ 2.0.*n2nt) .* dfv.n_H .* SA.abundance_for_element("Si")

# do Magnesium now
n1nt = SA.species_fraction.(dfv.T, Pe2, "Mg", ion="First")
n2nt = SA.species_fraction.(dfv.T, Pe2, "Mg", ion="Second")
n_e_Mg = (n1nt .+ 2.0.*n2nt) .* dfv.n_H .* SA.abundance_for_element("Mg")

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
ax1.set_yscale("log")
ax1.set_ylim(0.1,1.1)
ax1.set_xlabel(L"{\rm Height\ (km)}")
ax1.set_ylabel(L"{\rm Contribution\ to\ } n_{\rm e}")
ax1.legend(ncol=2)
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
pg = (SA.mH * sum(filter(!isnan, dfa.A .* dfa.Weight)))
pg = 2.2701e-24

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
