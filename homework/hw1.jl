# use packages
using Pkg; Pkg.activate(".")
using StellarAtmospheres
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib
using CSV
using DataFrames

# calc blackbody curve in frequency space for a few temperatures
T = [1000.0, 2000.0, 5000.0, 10000.0, 20000.0]
ν = range(1e3, 1e17, length=5000) # in Hz
I = zeros(length(ν), length(T))
for i in 1:length(T)
    I[:,i] = Bν.(ν, T[i])
end

# plot it
mpl.rc("font", size=14)
fig = plt.figure("Planck Function")
ax1 = plt.axes()
for i in 1:size(I,2)
    p1 = plot(ν, I[:,i], ls="-", linewidth=2, c=plt.cm.Dark2(i))
end
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel(L"B_\nu \ ({\rm erg}\ {\rm s}^{-1}\ {\rm cm}^{-2}\ {\rm Hz}^{-1}\ {\rm ster}^{-1})")
ax1.set_ylim(1e-25, 1e-2)
ax1.set_xlim(1e3, 1e17)
ax1.grid(color="gray", ls=":")
ax1.legend([string(temp) * " K" for temp in T])
plt.tight_layout()
fig.savefig("/Users/michael/Desktop/ASTRO530/figures/hw1_planck_function.pdf")

# write peak temperatures to CSV
λpeak = λmax.(T)
νpeak = νmax.(T)
df = DataFrame(temp=T, wave=λpeak, freq=νpeak)
CSV.write("/Users/michael/Desktop/ASTRO530/data/wien.tsv", df, delim="\t")
