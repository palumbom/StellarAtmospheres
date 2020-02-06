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
b = range(0.1, 25.0, length=100)

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

# find best result
am1 = argmin(int1_err); println(int1[am1])

# now do the transformation
u = range(-10, 2.5, length=100)
ab = (-1e5, 1e5)
int1_log = similar(u)

# calculate the integrand
integrand1 = expint_log.(1, u)
integrand2 = expint_log.(2, u)
integrand3 = expint_log.(3, u)

# plot the integrand
fig = plt.figure("Integrand")
ax1 = fig.add_subplot(111)
img = ax1.plot(u, integrand1, "-", label=L"n=1")
img = ax1.plot(u, integrand2, "-.", label=L"n=2")
img = ax1.plot(u, integrand3, ":", label=L"n=3")
ax1.set_xlabel(L"u")
ax1.set_ylabel(L"\log(10)10^u E_n(10^u)")
ax1.set_xlim(minimum(u), maximum(u))
ax1.set_ylim(1e-10, 1e0)
ax1.set_yscale("log")
ax1.legend()
fig.savefig(outdir*"hw4_integrand.pdf")
plt.clf(); plt.close()

# do the new integral
for i in eachindex(u)
    int1_log[i] = trap_int(x -> expint_log(1, u[i]), ab, ntrap=100, logx=false)
end
