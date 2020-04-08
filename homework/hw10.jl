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

# get stuff at τ=1 surface from VALIIIC
ind = 49
temp = dfv.T[ind]
Pg = dfv.Pg_Ptot[ind] .* dfv.Ptot[ind]
Pe = SA.P_from_nkt.(dfv.n_e[ind], temp)
ξ = dfv.V[ind] * 1000.0 * 100 # convert km/s -> cm/s
nH = dfv.n_H[ind]
ρ = dfv.ρ[ind]

# get neutral and ground fraction
f_neutral = SA.neutral_fraction(temp, Pe, "Na")
f_ground = 2.0/SA.calc_partition(temp, "Na") # equation 1.18

# input parameters for sodium
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")

# make LineParams object instances for NaD lines
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

# calculate sodium opacities as function of lambda
waves = range(5888.0, 5898.0, length=5000)
stime = SA.calc_SE.(waves, temp)
κ_na1 = SA.calc_α(waves, temp, Pe, Pg, ξ, NaD1) .* f_ground .* f_neutral .* na_a .* nH .* stime ./ ρ
κ_na2 = SA.calc_α(waves, temp, Pe, Pg, ξ, NaD2) .* f_ground .* f_neutral .* na_a .* nH .* stime ./ ρ
κcont = SA.κ_tot.(waves, temp, Pe, Pg)
κ_tot = κ_na1 + κ_na2 + κcont

# plot it
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(waves, κ_na1, "r-.", label=L"{\rm Na\ I\ D}_1")
ax1.plot(waves, κ_na2, "r--", label=L"{\rm Na\ I\ D}_2")
ax1.plot(waves, κcont, "b:", label=L"{\rm Continuous}")
ax1.plot(waves, κ_tot, "k-", label=L"{\rm Total}")
ax1.set_xlabel(L"{\rm Wavelength}\ ({\rm \AA})")
ax1.set_ylabel(L"{\rm Opacity}\ ({\rm cm}^2\ {\rm g}^{-1})")
ax1.legend()
fig.savefig(outdir*"hw10_kappa.pdf")
plt.clf(); plt.close()

# plot it
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_yscale("log")
ax1.plot(waves, κ_na1, "r-.", label=L"{\rm Na\ I\ D}_1")
ax1.plot(waves, κ_na2, "r--", label=L"{\rm Na\ I\ D}_2")
ax1.plot(waves, κcont, "b:", label=L"{\rm Continuous}")
ax1.plot(waves, κ_tot, "k-", label=L"{\rm Total}")
ax1.set_xlabel(L"{\rm Wavelength}\ ({\rm \AA})")
ax1.set_ylabel(L"{\rm Opacity}\ ({\rm cm}^2\ {\rm g}^{-1})")
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=2, mode="expand", borderaxespad=0.)
fig.savefig(outdir*"hw10_log_kappa.pdf")
plt.clf(); plt.close()
