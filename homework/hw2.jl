# use packages
using Pkg; Pkg.activate(".")
using StellarAtmospheres
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib

# array of temperatures and wavenumbers
T = [3000.0, 7000.0, 10000.0]
ν̃ = range(0.01, 12.0, length=1000)

# evaluate the planck function
Iν̃ = zeros(length(ν̃), length(T))
for i in 1:length(T)
    Iν̃[:,i] = Bν̃.(ν̃.*1e4, T[i])
end

# set some matplotlib params for pretty plots
mpl.rc("font", size=14)
mpl.rc("lines", linewidth=2, linestyle="-")
mpl.rc("axes", grid=true)
mpl.rc("grid", linestyle=":")
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
plt.ioff()

# linlin plot
fig = plt.figure("Planck Function")
ax1 = fig.add_subplot(111)
for i in 1:size(Iν̃,2)
    p1 = plot(ν̃, Iν̃[:,i], c=plt.cm.Dark2(i))
end
ax1.set_xlim(0.0, 12.0)
ax1.ticklabel_format(axis="both", style="sci")
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%1.0e"))
ax1.set_xlabel("Wavenumber " * L"({\rm \mu m}^{-1})")
ax1.set_ylabel(L"B_\tilde{\nu} \ ({\rm erg}\ {\rm \mu m}\ {\rm s}^{-1}\ {\rm cm}^{-2}\ {\rm ster}^{-1})")
ax1.legend([string(temp) * " K" for temp in T])
plt.tight_layout()
fig.savefig(outdir * "hw2_planck.pdf")
plt.clf()

# semilog plot
fig = plt.figure("Semilog Planck Function")
ax2 = fig.add_subplot(111)
for i in 1:size(Iν̃,2)
    p1 = plot(ν̃, Iν̃[:,i], c=plt.cm.Dark2(i))
end
ax2.set_xlim(0.0, 12.0)
ax2.set_xlabel("Wavenumber " * L"({\rm \mu m}^{-1})")
ax2.set_ylabel(L"B_\tilde{\nu} \ ({\rm erg}\ {\rm \mu m}\ {\rm s}^{-1}\ {\rm cm}^{-2}\ {\rm ster}^{-1})")
ax2.set_yscale("log")
ax2.legend([string(temp) * " K" for temp in T])
plt.tight_layout()
fig.savefig(outdir * "hw2_semilog_planck.pdf")
plt.clf()

# loglog plot
fig = plt.figure("Loglog Planck Function")
ax3 = fig.add_subplot(111)
for i in 1:size(Iν̃,2)
    p1 = plot(ν̃, Iν̃[:,i], c=plt.cm.Dark2(i))
end
ax3.set_xlabel("Wavenumber " * L"({\rm \mu m}^{-1})")
ax3.set_ylabel(L"B_\tilde{\nu} \ ({\rm erg}\ {\rm \mu m}\ {\rm s}^{-1}\ {\rm cm}^{-2}\ {\rm ster}^{-1})")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(10^-1.0, 10^1.2)
ax3.legend([string(temp) * " K" for temp in T])
plt.tight_layout()
fig.savefig(outdir * "hw2_loglog_planck.pdf")
plt.clf()
