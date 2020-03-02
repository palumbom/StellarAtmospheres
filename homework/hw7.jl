using Pkg; Pkg.activate(".")
using Revise
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# calculate opacities
λs = range(3000.0, 20000.0, length=1000)

# do Gray Fig. 8.5b first
T = 5143.0
Pe = exp10(1.08)
κHbf = SA.κ_H_bf.(λs, T, Pe)./Pe
κHff = SA.κ_H_ff.(λs, T, Pe)./Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe)./Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe)./Pe
κT = SA.κ_tot.(λs, T, Pe)./Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label="Total")
ax1.plot(λs, κHbf, "b--", label="H B-F")
ax1.plot(λs, κHff, "b:", label="H F-F")
ax1.plot(λs, κHmbf, "r--", label="H- B-F")
ax1.plot(λs, κHmff, "r:", label="H- F-F")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 5143\ {\rm K}", (16000, 6.9e-26))
ax1.annotate(L"\log(P_e) = 1.08", (15500, 6.1e-26))
fig.savefig(outdir * "hw7_85a.pdf")
plt.clf(); plt.close()

# do Gray Fig. 8.5b first
T = 6429.0
Pe = exp10(1.77)
κHbf = SA.κ_H_bf.(λs, T, Pe)./Pe
κHff = SA.κ_H_ff.(λs, T, Pe)./Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe)./Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe)./Pe
κT = SA.κ_tot.(λs, T, Pe)./Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label="Total")
ax1.plot(λs, κHbf, "b--", label="H B-F")
ax1.plot(λs, κHff, "b:", label="H F-F")
ax1.plot(λs, κHmbf, "r--", label="H- B-F")
ax1.plot(λs, κHmff, "r:", label="H- F-F")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 5143\ {\rm K}", (16000, 2.8e-26))
ax1.annotate(L"\log(P_e) = 1.77", (15500, 2.52e-26))
fig.savefig(outdir * "hw7_85b.pdf")
plt.clf(); plt.close()

# do Gray Fig. 8.5c now
T = 7715.0
Pe = exp10(2.50)
κHbf = SA.κ_H_bf.(λs, T, Pe)./Pe
κHff = SA.κ_H_ff.(λs, T, Pe)./Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe)./Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe)./Pe
κT = SA.κ_tot.(λs, T, Pe)./Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label="Total")
ax1.plot(λs, κHbf, "b--", label="H B-F")
ax1.plot(λs, κHff, "b:", label="H F-F")
ax1.plot(λs, κHmbf, "r--", label="H- B-F")
ax1.plot(λs, κHmff, "r:", label="H- F-F")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 7715\ {\rm K}", (16000, 3.5e-26))
ax1.annotate(L"\log(P_e) = 2.50", (15500, 3.1e-26))
fig.savefig(outdir * "hw7_85c.pdf")
plt.clf(); plt.close()

# do Gray Fig. 8.5d now
T = 11572.0
Pe = exp10(2.76)
κHbf = SA.κ_H_bf.(λs, T, Pe)./Pe
κHff = SA.κ_H_ff.(λs, T, Pe)./Pe
κHmbf = SA.κ_H_minus_bf.(λs, T, Pe)./Pe
κHmff = SA.κ_H_minus_ff.(λs, T, Pe)./Pe
κT = SA.κ_tot.(λs, T, Pe)./Pe

# plot it
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(λs, κT, "k-", label="Total")
ax1.plot(λs, κHbf, "b--", label="H B-F")
ax1.plot(λs, κHff, "b:", label="H F-F")
ax1.plot(λs, κHmbf, "r--", label="H- B-F")
ax1.plot(λs, κHmff, "r:", label="H- F-F")
ax1.set_xlabel(L"{\rm Wavelength}\ (\AA)")
ax1.set_ylabel(L"\kappa_\nu/P_e\ ({\rm cm}^2{\rm /H\ atom\ per\ dyne/cm}^2)")
ax1.legend(ncol=2, fontsize=10)
ax1.set_xticks(range(5000, 20000, length=4))
ax1.set_xlim(2000, 21000)
ax1.annotate(L"T = 11572\ {\rm K}", (16000, 1.1e-25))
ax1.annotate(L"\log(P_e) = 2.76", (15500, 1.0e-25))
fig.savefig(outdir * "hw7_85d.pdf")
plt.clf(); plt.close()
