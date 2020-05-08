using Pkg; Pkg.activate(".")
using Dierckx
using BenchmarkTools
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# interpolate departure coefficients and source functions
npoints = 3000
depc_df = SA.interp_nlte(SA.dfd, npoints=npoints)
source_df = SA.interp_nlte(SA.dfs, npoints=npoints)

# actuall pull out the height, b, S arrays
h = depc_df.h .* 1000.0 .* 100.0 # convert km -> cm
b = depc_df.b
S = source_df.S

# convert units of S to per Hz from per nm
# the factor of 10 was revealed to me in a dream
S *= 10.0 * 1e9 * 3e8/(λ2ν(5890.0*1e-8))^2

# get interpolated VALIIIc
val_new = SA.interp_valIIIc(npoints=npoints)

# assign VALIIIc variables for convenience
temp = val_new.T
ne = val_new.n_e
nH = val_new.n_H
Pe = SA.P_from_nkt.(ne, temp)
Pg = val_new.Pg_Ptot .* val_new.Ptot
ξ = val_new.V .* 1000.0 .* 100.0 # convert km/s -> cm/s
ρ = val_new.ρ
τ_500 = val_new.τ_500

# make LineParams object instances for NaD lines
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

# get the total opacity using departure coefficients
λs = range(5888.0, 5898.0, step=0.05)
κ_na1 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD1, dep=b)
κ_na2 = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD2, dep=b)
κcont = SA.κ_tot(λs, temp, Pe, Pg)
κ_tot = (κcont .+ κ_na1 .+ κ_na2)

# get LTE opacs
κ_na1_lte = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD1)
κ_na2_lte = κ_line(λs, temp, Pe, Pg, ξ, nH, ρ, NaD2)
κ_tot_lte = (κcont .+ κ_na1_lte .+ κ_na2_lte)

# calculate τν
τ_mid = elav(τ_500)
τν = calc_τν(κ_tot, ρ, h)
τν_lte = calc_τν(κ_tot_lte, ρ, h)

# plot calc tau against wavelength
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_yscale("log")
img = ax1.contourf(λs, τ_mid, log10.(τν), levels=range(-7.5, 4.5, length=9))
cbr = fig.colorbar(img)
cbr.set_label(L"\log\tau_\nu")
ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
ax1.set_ylabel(L"\tau_{500}")
ax1.set_ylim(1e-7, 5.0)
fig.savefig(outdir * "hw13_tau_image.pdf")
plt.clf(); plt.close()

# explicit Milne integral
Tsun = 5777.0
νs = λ2ν.(λs .* 1e-8)
the_flux = similar(νs)
the_flux_lte = similar(νs)
for i in eachindex(νs)
    ys = elav(S) .* SA.expint.(2, τν[:,i])
    ys_lte = SA.Bν.(νs[i], elav(temp)) .* SA.expint.(2, τν_lte[:,i])
    the_flux[i] = SA.trap_int(τν[:,i], ys)
    the_flux_lte[i] = SA.trap_int(τν_lte[:,i], ys_lte)
end

# convert to flux
the_flux .*= 2π
the_flux_lte .*= 2π

# plot it
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(λs, the_flux./1e-5, "k-", label=L"{\rm NLTE}")
ax1.plot(λs, the_flux_lte./1e-5, "r", alpha=0.75, label=L"{\rm LTE}")
ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
ax1.set_ylabel(L"\mathcal{F}_\nu^+\ " * SA.int_units_string())
ax1.legend()
fig.savefig(outdir * "hw13_emergent_flux.pdf")
plt.clf(); plt.close()

# also plot zoom ins
ind_mid = SA.searchsortednearest(λs, 5893.0)
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(λs[1:ind_mid], the_flux[1:ind_mid]./1e-5, "k-", label=L"{\rm NLTE}")
ax1.plot(λs[1:ind_mid], the_flux_lte[1:ind_mid]./1e-5, "r", alpha=0.75, label=L"{\rm LTE}")
ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
ax1.set_ylabel(L"\mathcal{F}_\nu^+\ " * SA.int_units_string())
ax1.legend()
fig.savefig(outdir * "hw13_emergent_d2.pdf")
plt.clf(); plt.close()

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(λs[ind_mid:end], the_flux[ind_mid:end]./1e-5, "k-", label=L"{\rm NLTE}")
ax1.plot(λs[ind_mid:end], the_flux_lte[ind_mid:end]./1e-5, "r", alpha=0.75, label=L"{\rm LTE}")
ax1.set_xlabel(L"{\rm Wavelength\ (\AA)}")
ax1.set_ylabel(L"\mathcal{F}_\nu^+\ " * SA.int_units_string())
ax1.legend()
fig.savefig(outdir * "hw13_emergent_d1.pdf")
plt.clf(); plt.close()

# get the mean intensity
core_ind = SA.searchsortednearest(λs, 5890.0)
the_J_wing = zeros(length(h)-1)
the_J_core = zeros(length(h)-1)
for i in eachindex(the_J_wing)
    ys = elav(S) .* SA.expint.(1, abs.(τν[:,ind_mid] .- τν[i,ind_mid]))
    ys2 = elav(S) .* SA.expint.(1, abs.(τν[:,core_ind] .- τν[i,core_ind]))
    inds = .!isinf.(ys)
    inds2 = .!isinf.(ys2)
    the_J_wing[i] = SA.trap_int(τν[:,ind_mid][inds], ys[inds])
    the_J_core[i] = SA.trap_int(τν[:,core_ind][inds], ys2[inds])
end
the_J_wing .*= 0.5
the_J_core .*= 0.5

# now replicate the rutten plot
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_yscale("log")
ax1.plot(h./1e5, SA.Bν.(νs[ind_mid], temp), "k:", label=L"B_\nu")
ax1.plot(h./1e5, S, "k", label=L"S_\nu")
ax1.plot(elav(h)./1e5, the_J_wing, "k--", label=L"J_{\rm wing}")
ax1.plot(elav(h)./1e5, the_J_core, "k-.", label=L"J_{\rm core}")
ax1.set_xlim(-100, 2000)
ax1.set_ylim(1e-7, 2e-4)
ax1.set_xlabel(L"{\rm Height\ (km)}")
ax1.set_ylabel(L"{\rm Source\ function}\ " * SA.int_units_string())
ax1.legend()
fig.savefig(outdir * "hw13_rutten_replica.pdf")
plt.clf(); plt.close()
