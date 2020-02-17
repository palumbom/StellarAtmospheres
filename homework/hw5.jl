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
ν̃ = range(0.01, 25.0, length=500)
ν = ν̃2ν.(ν̃.*1e4)
T = 8700.0
τs = (1e-5, 25.0)

# calculate the emergent flux
F1, err = ℱν₀(SνPlanck, τs, Teff=T, ν=ν, ntrap=5000, err=true)
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
ax1.annotate(L"T_{\rm eff} = 8700\ {\rm K}", (0.1, 5e-4), xycoords="data")
ax2.plot(ν̃, F1 - F2, "k-")
ax2.set_xlim(-0.5, 12.5)
ax2.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax2.set_ylabel(L"\mathcal{F}_\nu^{\rm Int.}(0) - \mathcal{F}_\nu^{\rm EB}(0)")
fig.savefig(outdir * "hw5_F_nu_0_8700.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# plot the error
fig = plt.figure(figsize=(8,7))
ax1 = fig.add_subplot()
ax1.plot(ν̃, abs.(err), "k-")
ax1.set_xlim(-0.5, 12.5)
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"{\rm Approximate\ Int.\ Error}  \ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.set_yscale("log")
ax1.set_ylim(1e-16, 1e-6)
fig.savefig(outdir * "hw5_int_error.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# second derivative of Sν
τ = range(0.0, 2.0, length=51)
d2 = similar(ν)
for i in eachindex(ν)
    f = t -> SνPlanck(t, Teff=T, ν=[ν[i]])[1]
    d2[i] = deriv2(f, τ)[26]
end

# plot the second derivative
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(ν̃, d2, "k-")
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\frac{d^2S_\nu(\tau)}{d\tau^2}|_{\tau = 1}\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.set_xlim(-0.5, 12.5)
fig.savefig(outdir * "hw5_dSdT.pdf")
plt.clf(); plt.close()

# calculate & plot flux at values of τ
τ = [0.1, 1.0, 2.0, 5.0, 10.0, 20.0]
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in eachindex(τ)
    lab = L"\tau=\ " * string(τ[i])
    ax1.plot(ν̃, ℱντ(SνPlanck, τ[i], (0.0, Inf), Teff=T, ν=ν, quad=true), c=plt.cm.Dark2(i), label=lab)
end
ax1.set_xlim(-0.5,16.0)
ax1.set_xlabel(L"\tilde{\nu}\ (\mu{\rm m}^{-1})")
ax1.set_ylabel(L"\mathcal{F}_\nu (\tau)\ ({\rm erg/s/cm}^2{\rm /Hz/ster})")
ax1.legend(ncol=2)
fig.savefig(outdir * "hw5_F_nu_tau_8700.pdf", bbox_inches="tight")
plt.clf(); plt.close()

# integrate ℱντ over ν
τ = range(0.1, 10.0, length=300)
Fτs = similar(τ)
for i in eachindex(τ)
    Fτs[i] = ℱτ(SνPlanck, τ[i], τs, Teff=T, ν=ν, quad=true)
end

# get Stefan Boltzmann
Fsb = SB(T)

# now plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(τ, Fτs, "k-", label="Integration")
ax1.plot(τ, fill(Fsb, (length(τ))), "r:", label="Stefan-Boltzmann")
ax1.set_yscale("log")
ax1.set_xlim(-1.0, τ[end]+1.0)
ax1.set_xlabel(L"\tau")
ax1.set_ylabel(L"\mathcal{F}(\tau)\ ({\rm erg/s/cm}^2{\rm /ster})")
ax1.legend()
fig.savefig(outdir*"hw5_F_tau_8700.pdf", bbox_inches="tight")
plt.clf(); plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(τ, abs.(Fτs .- Fsb)./Fsb, "k-")
ax1.set_xlim(-1.0, τ[end]+1.0)
ax1.set_yscale("log")
ax1.set_xlabel(L"\tau")
ax1.set_ylabel(L"|\mathcal{F}(\tau) - \mathcal{F}_{\rm SB} (\tau)| /\mathcal{F}_{\rm SB} (\tau)")
fig.savefig(outdir*"hw5_sb.pdf", bbox_inches="tight")
plt.clf(); plt.close()
