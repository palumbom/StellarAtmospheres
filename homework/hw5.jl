using Pkg; Pkg.activate(".")
using Revise
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# parameters for flux calc
ν̃ = range(0.01, 15.0, length=500)
ν = ν̃2ν.(ν̃.*1e4)
T = 8700.0
τs = (1e-10, 50.0)

# calculate the emergent flux
F1 = ℱν₀(SνPlanck, τs, Teff=T, ν=ν, ntrap=500)
F2 = ℱνEB(SνPlanck, Teff=T, ν=ν)

# plot the emergent flux
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(ν̃, F1, "k-", label=L"2\pi \int_0^\infty S_\nu (t_\nu) E_2(t_\nu) dt_\nu")
ax1.plot(ν̃, F2, "r:", label="Eddington-Barbier")
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\mathcal{F}_\nu (0)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.set_xlim(-0.5, 12.5)
ax1.set_ylim(1e-7, 1e-3)
ax1.set_yscale("log")
ax1.legend()
ax1.annotate(L"T_{\rm eff} = 8777\ {\rm K}", (0.1, 5e-4), xycoords="data")
fig.savefig(outdir * "hw5_F_nu_0_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# calculate & plot flux at values of τ
τ = [0.1, 2.0, 5.0, 10.0]
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in eachindex(τ)
    lab = L"\tau=\ " * string(τ[i])
    ax1.plot(ν̃, ℱντ(SνPlanck, τ[i], τs, Teff=T, ν=ν, ntrap=2000), label=lab)
end
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\mathcal{F}_\nu (\tau)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
# ax1.annotate(L"T_{\rm eff} = 8777\ {\rm K}", (0.2, 4.2e-4), xycoords="data")
ax1.legend(ncol=2)
fig.savefig(outdir * "hw5_F_nu_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()

ℱτ(SνPlanck, 1.0, τs, Teff=T, ν=ν, ntrap=1000)

# integrate ℱντ over ν
Fτs = similar(τ)
for i in eachindex(τ)
    Fτs[i] = ℱτ(SνPlanck, τ[i], τs, Teff=T, ν=ν, ntrap=300)
end

# now plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(τ, Fτs)
ax1.set_xlabel(L"\tau")
ax1.set_ylabel(L"\mathcal{F}(\tau)\ ({\rm erg/s/cm}^2{\rm /ster})")
ax1.set_xlim(τ[1], τ[end])
fig.savefig(outdir*"hw5_F_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()
