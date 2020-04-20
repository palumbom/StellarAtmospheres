using Pkg; Pkg.activate(".")
using Dierckx
using BenchmarkTools
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# get VALIIIc on interpolated grid
val_old = SA.dfv
val_new = SA.interp_valIIIc(npoints=100)

# old quantities
temp_old = val_old.T
ne_old = val_old.n_e
nH_old = val_old.n_H
Pe_old = SA.P_from_nkt.(ne_old, temp_old)
Pg_old = val_old.Pg_Ptot .* val_old.Ptot
ξ_old = val_old.V .* 1000.0 .* 100.0 # convert km/s -> cm/s
ρ_old = val_old.ρ
h_old = val_old.h .* 1000.0 .* 100.0 # convert km -> cm
τ_500_old = val_old.τ_500

# assign VALIIIc variables for convenience
temp = val_new.T
ne = val_new.n_e
nH = val_new.n_H
Pe = SA.P_from_nkt.(ne, temp)
Pg = val_new.Pg_Ptot .* val_new.Ptot
ξ = val_new.V .* 1000.0 .* 100.0 # convert km/s -> cm/s
ρ = val_new.ρ
h = val_new.h #.* 1000.0 .* 100.0 # convert km -> cm
τ_500 = val_new.τ_500

# wavelength range generator + continuum opacities
λs = range(5886.0, 5900.0, length=300)
λs = 5890.0
κcont = SA.κ_tot.(λs, temp, Pe, Pg)

# make LineParams object instances for NaD lines
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

# calculate line opacity
κ_na1 = zeros(length(temp))
κ_na2 = zeros(length(temp))
for i in eachindex(κ_na1)
    # get alpha
    α1 = SA.calc_α(λs, temp[i], Pe[i], Pg[i], ξ[i], NaD1)
    α2 = SA.calc_α(λs, temp[i], Pe[i], Pg[i], ξ[i], NaD2)

    # other quantities
    f_neutral = SA.neutral_fraction(temp[i], Pe[i], "Na")
    f_ground = 2.0/exp10(SA.calc_partition(temp[i], "Na")) # equation 1.18
    na_a = SA.abundance_for_element("Na")
    stime = SA.calc_SE.(λs, temp[i])

    κ_na1[i] = α1 * f_ground * f_neutral * na_a * nH[i] * stime / ρ[i]
    κ_na2[i] = α2 * f_ground * f_neutral * na_a * nH[i] * stime / ρ[i]
end

# total opacity
κ_tot = (κcont .+ κ_na1 .+ κ_na2)#'

# # calculate τν
# τ_mid = elav(τ_500)

# τν = calc_τν(κ_tot, ρ, h)

# # plot calc tau against wavelength
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.set_yscale("log")
# img = ax1.contourf(λs, τ_mid, log10.(τν))
# cbr = fig.colorbar(img)
# cbr.set_label(L"\log\tau_\nu")
# ax1.set_xlabel(L"{\rm Wavelength\ \AA}")
# ax1.set_ylabel(L"\tau_{500}")
# ax1.set_ylim(1e-7, 5.0)
# fig.savefig(outdir * "hw11_tau_image.pdf")
# plt.clf(); plt.close()

# # explicit
# Tsun = 5777.0
# νs = λ2ν.(λs .* 1e-8)
# the_flux = similar(λs)
# for i in eachindex(λs)
#     delta_τ = diff(τν[:,i])
#     sfunc = SA.SνPlanck.(νs[i], τν[:,i], Teff=Tsun) .* SA.expint.(2, τν[:,i])
#     the_flux[i] = 2π * sum(elav(sfunc) .* delta_τ)
# end

# # delta_τ = diff(τν, dims=1)
# # sfunc = SA.SνPlanck.(νs[i], τν[:,i], Teff=Tsun) .* SA.expint.(2, τν[:,i])

# # plot it
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.plot(λs, the_flux)
# ax1.set_xlabel(L"{\rm Wavelength\ \AA}")
# ax1.set_ylabel(L"\mathcal{F}_\nu(0)\ {\rm cgs}")
# fig.savefig(outdir * "hw11_emergent_flux.pdf")
# plt.clf(); plt.close()
