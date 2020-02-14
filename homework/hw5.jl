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
ν̃ = range(0.01, 25.0, length=1000)
ν = ν̃2ν.(ν̃.*1e4)
T = 8700.0
τs = (1e-10, 100.0)

# calculate the emergent flux
F1 = ℱν₀(SνPlanck, τs, Teff=T, ν=ν, ntrap=5000)
F2 = ℱνEB(SνPlanck, Teff=T, ν=ν)

# plot the emergent flux
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 7))
ax1.plot(ν̃, F1, "k-", label=L"2\pi \int_0^\infty S_\nu (t_\nu) E_2(t_\nu) dt_\nu")
ax1.plot(ν̃, F2, "r:", label="Eddington-Barbier")
ax1.set_ylabel(L"\mathcal{F}_\nu (0)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.set_xlim(-0.5, 12.5)
ax1.set_ylim(1e-7, 1e-3)
ax1.set_yscale("log")
ax1.legend()
ax1.annotate(L"T_{\rm eff} = 8777\ {\rm K}", (0.1, 5e-4), xycoords="data")
ax2.plot(ν̃, F1 - F2, "k-")
ax2.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax2.set_ylabel(L"{\rm Residual}\ \mathcal{F}_\nu (0)")
fig.savefig(outdir * "hw5_F_nu_0_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# calculate & plot flux at values of τ
τ = [0.1, 1.0, 2.0, 5.0, 10.0]
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in eachindex(τ)
    lab = L"\tau=\ " * string(τ[i])
    ax1.plot(ν̃, ℱντ(SνPlanck, τ[i], τs, Teff=T, ν=ν, ntrap=5000), label=lab)
end
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\mathcal{F}_\nu (\tau)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.legend(ncol=2)
fig.savefig(outdir * "hw5_F_nu_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# integrate ℱντ over ν
τ = range(0.1, 10.1, length=50)
Fτs = similar(τ)
for i in eachindex(τ)
    Fτs[i] = ℱτ(SνPlanck, τ[i], τs, Teff=T, ν=ν, ntrap=5000)
end

# get Stefan Boltzmann
Fsb = SB(T)

# now plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(τ, Fτs, "k-", label="Integration")
ax1.plot(τ, fill(Fsb, (length(τ))), "r:", label="Stefan-Boltzmann")
ax1.set_xlabel(L"\tau")
ax1.set_ylabel(L"\mathcal{F}(\tau)\ ({\rm erg/s/cm}^2{\rm /ster})")
ax1.set_xlim(-1.0, τ[end]+1.0)
ax1.legend()
fig.savefig(outdir*"hw5_F_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()
