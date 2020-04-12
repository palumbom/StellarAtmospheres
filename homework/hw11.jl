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

# get VALIIIc stuff
dfv = SA.dfv
temp = dfv.T
ne = dfv.n_e
nH = dfv.n_H
Pe = SA.P_from_nkt.(ne, temp)
Pg = dfv.Pg_Ptot .* dfv.Ptot
ξ = dfv.V .* 1000.0 .* 100 # convert km/s -> cm/s
ρ = dfv.ρ
h = dfv.h .* 1000.0 .* 100 # convert km -> cm
τ_500 = dfv.τ_500

# wavelength range generator + continuum opacities
λs = range(5886.0, 5900.0, length=1000)
κcont = SA.κ_tot(λs, temp, Pe, Pg)

# make LineParams object instances for NaD lines
m = 3.817541e-23
na_a = SA.abundance_for_element("Na")
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

# sodium opacities
κ_na1 = κ_line(λs, temp[49], Pe[49], Pg[49], ξ[49], nH[49], ρ[49], NaD1)
κ_na2 = κ_line(λs, temp[49], Pe[49], Pg[49], ξ[49], nH[49], ρ[49], NaD2)

# total opacity
κ_tot = (κcont .+ κ_na1 .+ κ_na2)'

# calculate τν
τ_mid = elav(τ_500)
τν = calc_τν(κ_tot, ρ, h)

# plot calc tau against wavelength
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_yscale("log")
img = ax1.contourf(λs, τ_mid, log10.(τν))
cbr = fig.colorbar(img)
cbr.set_label(L"\log\tau_\nu")
ax1.set_xlabel(L"{\rm Wavelength\ \AA}")
ax1.set_ylabel(L"\tau_{500}")
ax1.set_ylim(1e-7, 5.0)
fig.savefig(outdir * "hw11_tau_image.pdf")
plt.clf(); plt.close()

# now do emergent flux
Tsun = 5777.0
spl = Spline2D(τ_mid, λs, τν)#, kx=3, ky=3)
τbounds = (minimum(τ_500), maximum(τ_500))

# extend spline method to handle tuples
(spl::Spline2D)(t::Tuple{T,T}) where T<:Real = spl(t[1], t[2])

# finer grid
λnew = λs
τnew = range(minimum(τ_500), maximum(τ_500), length=1000)
τνnew = map(spl, ((x,y) for x in τnew, y in λnew))

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.set_yscale("log")
img = ax1.contourf(λnew, τnew, log10.(τνnew))
cbr = fig.colorbar(img)
cbr.set_label(L"\log\tau_\nu")
ax1.set_xlabel(L"{\rm Wavelength\ \AA}")
ax1.set_ylabel(L"\tau_{500}")
ax1.set_ylim(1e-7, 5.0)
fig.savefig(outdir * "hw11_tau_image2.pdf")
plt.clf(); plt.close()

the_flux = similar(λs)
for i in eachindex(λs)
    the_flux[i] = ℱν₀_line(λ2ν(λs[i]*1e-8), spl, τbounds, Teff=Tsun)
end

plt.plot(λs, the_flux); plt.show()
