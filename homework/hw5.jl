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
ν̃ = range(0.01, 12.0, length=1000)
ν = ν̃2ν.(ν̃.*1e4)
T = 8700.0
τs = (1e-10, 100.0)

# calculate the emergent flux
F8777 = ℱν₀(SνPlanck, τs, Teff=T, ν=ν, ntrap=500)

# plot the emergent flux
ν̃lim = [0.0, 12.0]
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(ν̃, F8777, label="Integral")
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\mathcal{F}_\nu (0)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.set_xlim(ν̃lim)
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
ax1.set_xlim(ν̃lim)
ax1.set_ylim(-0.0005, 0.0005)
ax1.annotate(L"T_{\rm eff} = 8777\ {\rm K}", (0.2, 4.2e-4), xycoords="data")
ax1.legend(ncol=2)
fig.savefig(outdir * "hw5_F_nu_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()

derp = ℱντ(SνPlanck, 5.0, τs, Teff=T, ν=[3.53e13], ntrap=500)

# integrate ℱντ over ν
τ = range(0.01, 10.01, length=25)
Fτs = similar(τ)
for i in eachindex(τ)
    Fτs[i] = ℱτ(SνPlanck, τ[i], τs, Teff=T, ν=ν, ntrap=500)
end

# function fancy_label(num::T) where T<:Real
#     s = split(string(num), "e")
#     n = string(s[1])
#     e = string(s[2])
#     return L"$(@sprintf("%.2f", n)) \times $(@sprintf("%.2",e))"
# end

# now plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(τ, Fτs)
ax1.set_xlabel(L"\tau")
ax1.set_ylabel(L"\mathcal{F}(\tau)\ ({\rm erg/s/cm}^2{\rm /ster})")
ax1.set_xlim(τ[1], τ[end])
# yts = ax1.get_yticks()
# ax1.set_yticks(yts)
# ax1.set_yticklabels(yts)
# ax1.ticklabel_format(axis="y", style="sci", useOffset=false)
# # ax1.ticklabel_format(axis="y", style="sci", useOffset=false)
# ax1.set_yticklabels(ax1.get_yticks())
fig.savefig(outdir*"hw5_F_tau_8777.pdf", bbox_inches="tight")
plt.clf(); plt.close()
