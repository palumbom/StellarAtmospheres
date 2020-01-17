# use packages
using Pkg; Pkg.activate(".")
using StellarAtmospheres
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib

# array of temperatures and wavenumbers
T = [3000.0, 7000.0, 10000.0]
λ = range(0.01, 12.0, length=1000)   # !~ in microns ~!
ν̃ = λ2ν̃.(λ)

# evaluate the planck function
Iν̃ = zeros(length(λ), length(T))
for i in 1:length(T)
    Iν̃[:,i] = Bν̃.(ν̃.*1e4, T[i])
end

# make some plots
mpl.rc("font", size=14)
mpl.rc("lines", linewidth=2, linestyle="-")
mpl.rc("axes", grid=true
fig = plt.figure("Planck Function", figsize=(4,12))
ax1 = fig.add_subplot(3, 1, 1)
for i in 1:size(Iν̃,2)
    p1 = plot(ν̃, Iν̃[:,i], c=plt.cm.Dark2(i))
end
ax1.legend([string(temp) * " K" for temp in T])

ax2 = fig.add_subplot(3, 1, 2)
for i in 1:size(Iν̃,2)
    p1 = plot(ν̃, Iν̃[:,i], c=plt.cm.Dark2(i))
end
ax2.set_yscale("log")

ax3 = fig.add_subplot(3, 1, 3)
for i in 1:size(Iν̃,2)
    p1 = plot(ν̃, Iν̃[:,i], c=plt.cm.Dark2(i))
end
ax3.set_xscale("log")
ax3.set_yscale("log")
