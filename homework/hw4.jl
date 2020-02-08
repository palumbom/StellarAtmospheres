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
    int1[c] = trap_int(x -> expint(1, x), (a,b[j]), ntrap=ntrap[i], logx=true)
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
println(@benchmark trap_int(x -> expint(1, x), ab, ntrap=ntrap, logx=true))

# compare Hν(O) calculated by different means
τs = (1e-10, 100.0)
a01 = [1.0, 1.0]
a2 = range(-0.1, 0.1, length=100)

# compute via expint
Ha = similar(a2)
for i in eachindex(a2)
    Ha[i] = Hν₀(τs, a01..., a2[i], ntrap=1000)
end

# compute via integration Iν(0,μ)*μ dμ
Hb = similar(a2)
for i in eachindex(a2)
    Hb[i] = Hν₀(a01..., a2[i], ntrap=1000)
end

# compute via EB approx for Hν(0)
Hd = similar(a2)
for i in eachindex(a2)
    Hd[i] = HνEB(a01..., a2[i])
end

# now visualize the result
fig = plt.figure("Emergent Flux")
ax1 = fig.add_subplot(111)
ax1.plot(a2, Ha, "k-", label=L"\frac{1}{2} \int_0^\infty S_\nu (\tau_\nu) E_2(t_\nu) dt_\nu")
ax1.plot(a2, Hb, "k-.", label=L"\frac{1}{2} \int_0^1 I_\nu(0,\mu)\mu d\mu")
ax1.plot(a2, Hd, "k:", label=L"\frac{1}{2} S_\nu(\tau_\nu = \frac{2}{3})")
ax1.set_xlabel("Quadratic Coefficient " * L"a_2")
ax1.set_ylabel(L"H_\nu(0)")
ax1.set_xlim(a2[1], a2[end])
ax1.legend()
fig.savefig(outdir*"hw4_eddington_flux.pdf")
plt.clf(); plt.close()
