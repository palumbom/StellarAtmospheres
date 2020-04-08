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
nH = dfv.n_H
Pe = SA.P_from_nkt.(ne, temp)
Pg = dfv.Pg_Ptot .* dfv.Ptot
ξ = dfv.V .* 1000.0 .* 100 # convert km/s -> cm/s
ρ = dfv.ρ
h = (dfv.h .+ 75) .* 1000.0 .* 100 # convert km -> cm
τ_500 = dfv.τ_500

# get neutral and ground fraction
f_neutral = SA.neutral_fraction.(temp, Pe, "Na")
f_ground = 2.0 ./ SA.calc_partition.(temp, "Na") # equation 1.18

# wavelength range generator + continuum opacities
λs = range(5888.0, 5898.0, length=500)
κcont = SA.κ_tot(λs, temp, Pe, Pg)

# make LineParams object instances for NaD lines
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

# sodium opacities
κ_na1 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD1)
κ_na2 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD2)

# total opacity
κ_tot = κcont .+ κ_na1 .+ κ_na2

# get Δh and midpoints for rho and tau
Δh = -diff(h)
ρ_mid = elav(ρ)
τ_mid = elav(τ_500)

# get dtau at each level
dτ = zeros(length(λs), length(Δh))
for i in eachindex(Δh)
    dτ[:, i] = κ_tot[:, i] .* ρ_mid[i]' .* Δh[i]'
end

# cumulative sum (remember to transpose)
τν = cumsum(dτ', dims=1)

# plot calc tau against wavelength
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_yscale("log")
img = ax1.contourf(λs, τ_mid, log10.(τν))
cbr = fig.colorbar(img)
cbr.set_label(L"\log\tau_\nu")
ax1.set_xlabel(L"{\rm Wavelength\ \AA}")
ax1.set_ylabel(L"\tau_{500}")
ax1.set_ylim(1e-5, 5.0)
fig.savefig(outdir * "hw11_tau_image.pdf")
plt.clf(); plt.close()

# get NaD lines
# flux = ℱν₀()
