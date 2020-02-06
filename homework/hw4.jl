using Pkg; Pkg.activate(".")
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib;
import Bridge.expint

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# parameters for integrating expn function
ntrap = range(10, 10000, step=100)
a = 1e-10
b = range(1e1, 1e4, length=100)

# allocate memory
int1 = zeros(length(ntrap), length(b))
int2 = zeros(length(ntrap), length(b))
int3 = zeros(length(ntrap), length(b))

# do the integration
for c in CartesianIndices(int1)
    i,j = Tuple(c)
    int1[c] = trap_int(x -> expint(1, x), (a,b[j]), ntrap=ntrap[i], logx=true)
    int2[c] = trap_int(x -> expint(2, x), (a,b[j]), ntrap=ntrap[i], logx=true)
    int3[c] = trap_int(x -> expint(3, x), (a,b[j]), ntrap=ntrap[i], logx=true)
end

int1_err = log10.(abs.(int1 .- 1.0) ./ 1.0)
int2_err = log10.(abs.(int1 .- 0.5) ./ 0.5)
int3_err = log10.(abs.(int1 .- (1.0/3.0)) ./ (1.0/3.0))

# plot it
extent = [1e1, 1e4, 10, 1000]
fig = plt.figure("Integration Error")
ax1 = fig.add_subplot(111)
img = ax1.imshow(int1_err, cmap="plasma", extent=extent,
                 interpolation="nearest", origin="lower")
cbr = fig.colorbar(img)
cbr.set_label("log(Relative Integration Error)")
ax1.set_xlabel("Upper Bound of Integration")
ax1.set_ylabel("Number of Trapezoids")
ax1.set_aspect("auto")
fig.savefig(outdir*"hw4_integration_error.pdf")
plt.clf(); plt.close()
