# use packages
using Pkg; Pkg.activate(".")
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
using Conda; # Conda.add("astropy")
using PyCall; apm = pyimport("astropy.modeling.models"); u = pyimport("astropy.units")
using PyPlot; plt = PyPlot; mpl = matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# array of temperatures and wavenumbers
T = [3000.0, 7000.0, 10000.0]
ν̃ = range(0.01, 20.0, length=5000)

# evaluate the planck function
Iν = zeros(length(ν̃), length(T))
for i in 1:length(T)
    Iν[:,i] = Bν.(ν̃2ν.(ν̃.*1e4), T[i])
end

# make labels for legend
Tlabels = [string(temp) * " K" for temp in T]

# linlin plot
fig = plt.figure("Planck Function")
ax1 = fig.add_subplot(111)
for i in 1:size(Iν,2)
    p1 = plot(ν̃, Iν[:,i], c=plt.cm.Dark2(i))
end
ax1.set_xlim(0.0, 12.0)
ax1.ticklabel_format(axis="both", style="sci")
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%1.0e"))
ax1.set_xlabel(SA.ν̃_string())
ax1.set_ylabel(SA.Bν_string())
ax1.legend(Tlabels)
fig.savefig(outdir * "hw2_planck.pdf")
plt.clf()

# semilog plot
fig = plt.figure("Semilog Planck Function")
ax2 = fig.add_subplot(111)
for i in 1:size(Iν,2)
    p1 = plot(ν̃, Iν[:,i], c=plt.cm.Dark2(i))
end
ax2.set_xlabel(SA.ν̃_string())
ax2.set_ylabel(SA.Bν_string())
ax2.set_yscale("log")
ax2.set_xlim(0.0, 12.0)
ax2.set_ylim(1e-15, 1e-3)
ax2.legend(Tlabels)
fig.savefig(outdir * "hw2_semilog_planck.pdf")
plt.clf()

# loglog plot
fig = plt.figure("Loglog Planck Function")
ax3 = fig.add_subplot(111)
for i in 1:size(Iν,2)
    p1 = plot(ν̃, Iν[:,i], c=plt.cm.Dark2(i))
end
ax3.set_xlabel(SA.ν̃_string())
ax3.set_ylabel(SA.Bν_string())
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(10^-1.0, 10^1.2)
ax3.set_ylim(1e-15, 1e-3)
ax3.legend(Tlabels)
fig.savefig(outdir * "hw2_loglog_planck.pdf")
plt.clf()

# compare to astropy
PyBB = apm.BlackBody(temperature=3000.0 * u.K)(ν̃2ν.(ν̃.*1e4))

# analytical integral
int_true = (2*SA.h/SA.c^3) * (SA.kB/SA.h)^4 * (π^4/15) * 7500.0^4 * 1e-4

# integrate the function simply at first
x = range(0.01, 100.0, length=10000)
y = Bν.(ν̃2ν.(x.*1e4), 7500.0)
int_num1 = trap_int(x, y)

# now use a more complex implementation
# we need to make an anonymous function with one input
f = z -> Bν(ν̃2ν.(z.*1e4), 7500.0)
int_num2, err = trap_int(f, (0.01, 100.0), 10000, err=true)
