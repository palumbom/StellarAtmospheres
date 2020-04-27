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
val_new = SA.interp_valIIIc(npoints=1000)

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
h = val_new.h .* 1000.0 .* 100.0 # convert km -> cm
τ_500 = val_new.τ_500

# wavelength range generator + continuum opacities
λs = range(5888.0, 5898.0, step=0.1)
κcont = SA.κ_tot(λs, temp, Pe, Pg)

# make LineParams object instances for NaD lines
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

κ_na1 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD1)
κ_na2 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD2)

# total opacity
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
fig.savefig(outdir * "hw11_tau_image.pdf")
plt.clf(); plt.close()

# explicit Milne integral
Tsun = 5777.0
νs = λ2ν.(λs .* 1e-8)
the_flux = similar(νs)
for i in eachindex(νs)
    ys = SA.SνPlanck.(νs[i], τν[:,i], Teff=Tsun) .* SA.expint.(2, τν[:,i])
    the_flux[i] = SA.trap_int(τν[:,i], ys)
end

# plot it
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(λs, the_flux)
ax1.set_xlabel(L"{\rm Wavelength\ \AA}")
ax1.set_ylabel(L"\mathcal{F}_\nu(0)\ {\rm cgs}")
fig.savefig(outdir * "hw11_emergent_flux.pdf")
plt.clf(); plt.close()


# plot stuff to diagnose
ind = findfirst(λs .== 5890.0)
the_κ = κ_tot[:,ind]
the_τ = τν[:,ind]
the_S = SA.SνPlanck.(νs[ind], the_τ, Teff=Tsun)
the_e = SA.expint.(2, the_τ)
the_i = the_S .* the_e

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(h, the_κ)
ax1.set_xlabel("height (cm)")
ax1.set_ylabel("total opacity")
fig.savefig("/Users/michael/Desktop/opac.pdf")
plt.clf(); plt.close()

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(elav(h), the_τ)
ax1.set_xlabel("height midpoint (cm)")
ax1.set_ylabel("tau nu")
fig.savefig("/Users/michael/Desktop/tau.pdf")
plt.clf(); plt.close()

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(elav(h), the_S)
ax1.set_xlabel("height midpoint (cm)")
ax1.set_ylabel("source func")
fig.savefig("/Users/michael/Desktop/source.pdf")
plt.clf(); plt.close()

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(elav(h), the_e)
ax1.set_xlabel("height midpoint (cm)")
ax1.set_ylabel("expint")
fig.savefig("/Users/michael/Desktop/expint.pdf")
plt.clf(); plt.close()

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(elav(h), the_i)
ax1.set_xlabel("height midpoint (cm)")
ax1.set_ylabel("milne integrand")
fig.savefig("/Users/michael/Desktop/integrand.pdf")
plt.clf(); plt.close()


# confusogram
ind2 = findfirst(isapprox.(τ_500, 1.0, atol=1e-1))
fig, axs = plt.subplots(nrows=2, ncols=2)
axs[1,1].plot(νs, κ_tot[ind2,:])
axs[1,1].set_xlabel("frequency")
axs[1,1].set_ylabel("total opac")
axs[1,2].plot(elav(h), τν[:,ind])
axs[1,2].set_xlabel("height")
axs[1,2].set_ylabel("tau nu")
axs[2,1].plot(νs, the_flux)
axs[2,1].set_xlabel("frequency")
axs[2,1].set_ylabel("emergent flux")
axs[2,2].plot(elav(h), the_S)
axs[2,2].set_xlabel("height")
axs[2,2].set_ylabel("source function")
fig.savefig("/Users/michael/Desktop/confuse.pdf")
plt.clf(); plt.close()
