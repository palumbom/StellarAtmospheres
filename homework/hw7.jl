using Pkg; Pkg.activate(".")
using Revise
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()


# Phi as function of T
Ts = range(1000.0, 12000.0, length=1000)
Φ1 = SA.ΦT.(Ts, "Fe", ion=:First)
Φ2 = SA.ΦT.(Ts, "Fe+", ion=:First)

# adjust to fraction
F0 = (1.0./Φ1)./((1.0./Φ1) .+ 1.0 .+ Φ2)
F1 = 1.0./((1.0./Φ1) .+ 1.0 .+ Φ2)
F2 = Φ2./((1.0./Φ1) .+ 1.0 .+ Φ2)

# now plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(Ts, F0, "k--", label=L"{\rm FeI}")
ax1.plot(Ts, F1, "k-", label=L"{\rm FeII}")
ax1.plot(Ts, F2, "k-.", label=L"{\rm FeIII}")
ax1.set_xlabel(L"{\rm Temperature\ (K)}")
ax1.set_ylabel(L"{\rm Fraction\ in\ ionization\ stage}")
ax1.legend()
ax1.set_xlim(1000,11000)
fig.savefig(outdir * "hw7_1.10.pdf")
plt.clf(); plt.close()

# calculate opacities
λs = range(3000.0, 20000.0, length=1000)

# do Gray Fig. 8.5a first
T = 5143.0
Pe = exp10(1.08)
SE = (1.0 .- exp10.(-1.2398e4./λs.*SA.temp_to_theta(T)))
PT = (1.0 + SA.ΦT(T, "H")/Pe)
κHbf = SA.κ_H_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHff = SA.κ_H_ff.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe) ./ PT .* 1e26 ./ Pe
κT = SA.κ_tot.(λs, T, Pe) .* 1e26 ./ Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label=L"{\rm Total}")
ax1.plot(λs, κHbf, "b--", label=L"{\rm HI\ BF}")
ax1.plot(λs, κHff, "b:", label=L"{\rm HI\ FF}")
ax1.plot(λs, κHmbf, "r--", label=L"{\rm H}^-\ {\rm BF}")
ax1.plot(λs, κHmff, "r:", label=L"{\rm H}^-\ {\rm FF}")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e \times 10^{-26}\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 5143\ {\rm K}", (16000, 6.9))
ax1.annotate(L"\log(P_e) = 1.08", (15500, 6.1))
fig.savefig(outdir * "hw7_85a.pdf")
plt.clf(); plt.close()

# do Gray Fig. 8.5b now
T = 6429.0
Pe = exp10(1.77)
SE = (1.0 .- exp10.(-1.2398e4./λs.*SA.temp_to_theta(T)))
PT = (1.0 + SA.ΦT(T, "H")/Pe)
κHbf = SA.κ_H_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHff = SA.κ_H_ff.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe) ./ PT .* 1e26 ./ Pe
κT = SA.κ_tot.(λs, T, Pe) .* 1e26 ./ Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label=L"{\rm Total}")
ax1.plot(λs, κHbf, "b--", label=L"{\rm HI\ BF}")
ax1.plot(λs, κHff, "b:", label=L"{\rm HI\ FF}")
ax1.plot(λs, κHmbf, "r--", label=L"{\rm H}^-\ {\rm BF}")
ax1.plot(λs, κHmff, "r:", label=L"{\rm H}^-\ {\rm FF}")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e \times 10^{-26}\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 5143\ {\rm K}", (16000, 2.8))
ax1.annotate(L"\log(P_e) = 1.77", (15500, 2.52))
fig.savefig(outdir * "hw7_85b.pdf")
plt.clf(); plt.close()

# do Gray Fig. 8.5c now
T = 7715.0
Pe = exp10(2.50)
SE = (1.0 .- exp10.(-1.2398e4./λs.*SA.temp_to_theta(T)))
PT = (1.0 + SA.ΦT(T, "H")/Pe)
κHbf = SA.κ_H_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHff = SA.κ_H_ff.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe) ./ PT .* 1e26 ./ Pe
κT = SA.κ_tot.(λs, T, Pe) .* 1e26 ./ Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label=L"{\rm Total}")
ax1.plot(λs, κHbf, "b--", label=L"{\rm HI\ BF}")
ax1.plot(λs, κHff, "b:", label=L"{\rm HI\ FF}")
ax1.plot(λs, κHmbf, "r--", label=L"{\rm H}^-\ {\rm BF}")
ax1.plot(λs, κHmff, "r:", label=L"{\rm H}^-\ {\rm FF}")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e \times 10^{-26}\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 7715\ {\rm K}", (16000, 3.5))
ax1.annotate(L"\log(P_e) = 2.50", (15500, 3.1))
fig.savefig(outdir * "hw7_85c.pdf")
plt.clf(); plt.close()

# do Gray Fig. 8.5d now
T = 11572.0
Pe = exp10(2.76)
SE = (1.0 .- exp10.(-1.2398e4./λs.*SA.temp_to_theta(T)))
PT = (1.0 + SA.ΦT(T, "H")/Pe)
κHbf = SA.κ_H_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHff = SA.κ_H_ff.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe) .* SE ./ PT .* 1e26 ./ Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe) ./ PT .* 1e26 ./ Pe
κT = SA.κ_tot.(λs, T, Pe) .* 1e26 ./ Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label=L"{\rm Total}")
ax1.plot(λs, κHbf, "b--", label=L"{\rm HI\ BF}")
ax1.plot(λs, κHff, "b:", label=L"{\rm HI\ FF}")
ax1.plot(λs, κHmbf, "r--", label=L"{\rm H}^-\ {\rm BF}")
ax1.plot(λs, κHmff, "r:", label=L"{\rm H}^-\ {\rm FF}")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e \times 10^{-26}\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 11572\ {\rm K}", (3000, 6.10))
ax1.annotate(L"\log(P_e) = 2.76", (2700, 5.40))
fig.savefig(outdir * "hw7_85d.pdf")
plt.clf(); plt.close()
