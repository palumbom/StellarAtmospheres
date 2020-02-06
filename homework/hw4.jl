using Pkg; Pkg.activate(".")
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib;

# pycall the exponential function
using Conda; ENV["CONDA_JL_HOME"] = "/usr/local/anaconda3/envs/conda_jl";
using PyCall; expn = pyimport("scipy.special").expn;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# parameters for integrating expn function
ntrap = 1000
ab = (1e-5, 1e5)

# do the integration
int1 = trap_int(x -> expn(1, x), ab, ntrap=ntrap, logx=true)
int2 = trap_int(x -> expn(2, x), ab, ntrap=ntrap, logx=true)
int3 = trap_int(x -> expn(3, x), ab, ntrap=ntrap, logx=true)
