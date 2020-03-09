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
