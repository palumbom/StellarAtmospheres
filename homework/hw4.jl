using Pkg; Pkg.activate(".")
using Revise
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib;
import Bridge.expint

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# parameters for integrating expn function
ntrap = range(2, 100, step=1)
a = 1e-10
b = range(0.1, 1e2, length=5000)

# allocate memory
int1 = zeros(length(ntrap), length(b))

# do the integration
for c in CartesianIndices(int1)
    i,j = Tuple(c)
    int1[c] = trap_int(x -> x*expint(1, x), (a,b[j]), ntrap=ntrap[i], logx=true)
end

int1_err = log10.(abs.(int1 .- 1.0) ./ 1.0)

# plot it
extent = [b[1], b[end], ntrap[1], ntrap[end]]
fig = plt.figure("Integration Error")
ax1 = fig.add_subplot(111)
img = ax1.imshow(int1_err, cmap="plasma", extent=extent,
                 interpolation="nearest", origin="lower")
cbr = fig.colorbar(img)
cbr.set_label("log(Relative Integration Error)")
ax1.set_xlabel("Upper Bound of Integration")
ax1.set_ylabel("Number of Trapezoids")
ax1.set_aspect("auto")
fig.savefig(outdir*"hw4_expn1_error.pdf")
plt.clf(); plt.close()

# best int
am = argmin(int1_err)
ab = (1e-10, b[am[2]])
ntrap = ntrap[am[1]]
@benchmark trap_int(x -> x*expint(1, x), ab, ntrap=ntrap, logx=true)
