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
val_new = SA.interp_valIIIc(npoints=3200)

# assign VALIIIc variables for convenience
temp = val_new.T
ne = val_new.n_e
nH = val_new.n_H
Pe = SA.P_from_nkt.(ne, temp)
Pg = val_new.Pg_Ptot .* val_new.Ptot
ξ = val_new.V .* 1000.0 .* 100.0 # convert km/s -> cm/s
ρ = val_new.ρ
h = val_new.h .* 1000.0 .* 100.0 # convert km -> cm
τ_500 = val_new.τ_500

# make LineParams object instances for NaD lines
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

# get the total opacity
λs = range(5888.0, 5898.0, step=0.1)
κ_na1 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD1)
κ_na2 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD2)
κcont = SA.κ_tot(λs, temp, Pe, Pg)
κ_tot = (κcont .+ κ_na1 .+ κ_na2)

# calculate τν
τ_mid = elav(τ_500)
τν = calc_τν(κ_tot, ρ, h)

# plot calc tau against wavelength
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_yscale("log")
img = ax1.contourf(λs, τ_mid, log10.(τν), levels=range(-7.5, 4.5, length=9))
cbr = fig.colorbar(img)
cbr.set_label(L"\log\tau_\nu")
ax1.set_xlabel(L"{\rm Wavelength\ \AA}")
ax1.set_ylabel(L"\tau_{500}")
ax1.set_ylim(1e-7, 5.0)
# fig.savefig(outdir * "hw11_tau_image.pdf")
plt.clf(); plt.close()

# explicit Milne integral
Tsun = 5777.0
νs = λ2ν.(λs .* 1e-8)
the_flux = similar(νs)
for i in eachindex(νs)
    ys = SA.Bν.(νs[i], elav(temp)) .* SA.expint.(2, τν[:,i])
    the_flux[i] = SA.trap_int(τν[:,i], ys)
end
the_flux .*= 2π

# plot it
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(λs, the_flux./1e-5)
ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
ax1.set_ylabel(L"\mathcal{F}_\nu^+\ (10^{-5}\ {\rm ergs\ s}^{-1}\ {\rm cm}^{-2}\ {\rm Hz}^{-1}) ")
# fig.savefig(outdir * "hw11_emergent_flux.pdf")
# plt.clf(); plt.close()
