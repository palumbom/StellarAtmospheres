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
ν̃ = range(0.01, 12.0, length=5000)
ν = ν̃2ν.(ν̃.*1e4)
T = 8700.0
τs = (1e-10, 50.0)

# calculate the emergent flux
F8777 = ℱν₀(SνPlanck, τs, Teff=T, ν=ν, ntrap=500)

# plot the emergent flux
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(ν̃, F8777, label="Integral")
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\mathcal{F}_\nu (0)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.set_ylim(1e-7, 1e-3)
ax1.set_yscale("log")
ax1.legend()
ax1.annotate(L"T_{\rm eff} = 8777\ {\rm K}", (0.1, 5e-4), xycoords="data")
fig.savefig(outdir * "hw5_F_nu_0_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# calculate & plot flux at values of τ
τ = [0.01, 0.1, 0.5, 0.67, 1.0, 2.0, 3.0, 10.0]
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in eachindex(τ)
    lab = L"\tau=\ " * string(τ[i])
    ax1.plot(ν̃, ℱντ(SνPlanck, τ[i], τs, Teff=T, ν=ν, ntrap=500), label=lab)
end
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\mathcal{F}_\nu (\tau)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.set_xlim(ν̃[1], ν̃[end])
ax1.set_ylim(-0.0005, 0.0005)
ax1.annotate(L"T_{\rm eff} = 8777\ {\rm K}", (0.2, 4.2e-4), xycoords="data")
ax1.legend(ncol=2)
fig.savefig(outdir * "hw5_F_nu_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# integrate ℱντ over ν
τ = range(0.01, 10.01, length=25)
Fτs = similar(τ)
for i in eachindex(τ)
    Fτs[i] = ℱτ(SνPlanck, τ[i], τs, Teff=T, ν=ν, ntrap=500)
end

# now plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(τ, Fτs)
ax1.set_xlabel(L"\tau")
ax1.set_ylabel(L"\mathcal{F}(\tau)\ ({\rm erg/s/cm}^2{\rm /ster})")
ax1.set_xlim(τ[1], τ[end])
ax1.ticklabel_format(style="plain")
fig.savefig(outdir*"hw5_F_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()
