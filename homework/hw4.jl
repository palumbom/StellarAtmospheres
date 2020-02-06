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
ntrap = range(100, 10000, step=100)
a = 1e-10
b = range(1e1, 1e10, length=100)

# allocate memory
int1 = zeros(100, 100)
int2 = zeros(100, 100)
int3 = zeros(100, 100)

# do the integration
for c in CartesianIndices(int1)
    i,j = Tuple(c)
    int1[c] = trap_int(x -> expint(1, x), (a,b[i]), ntrap=ntrap[j], logx=true))
    int2[c] = trap_int(x -> expint(2, x), (a,b[i]), ntrap=ntrap[j], logx=true))
    int3[c] = trap_int(x -> expint(3, x), (a,b[i]), ntrap=ntrap[j], logx=true))
end

