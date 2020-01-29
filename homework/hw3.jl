using Pkg; Pkg.activate(".")
using Revise
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
using PyPlot; plt = PyPlot; mpl = matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# evaulate the source functions
τ  = range(0.0, 10.0, length=1000)
an = [1.0, 1.0, 1.0]
S1 = Sν_lin.(τ, an...)
S2 = Sν_quad.(τ, an...)

# plot the source functions
fig = plt.figure("Source Functions")
ax1 = fig.add_subplot(111)
ax1.plot(τ, S1, label="Linear")
ax1.plot(τ, S2, label="Quadratic")
ax1.set_xlabel(L"\tau_\nu")
ax1.set_ylabel(L"S_\nu(\tau_\nu)")
ax1.set_xlim(0,10)
ax1.legend()
fig.savefig(outdir*"hw3_source.pdf")
plt.clf(); plt.close()

# evaluate the contribution function
C1 = Cν.(τ, Sν_lin, an...)
C2 = Cν.(τ, Sν_quad, an...)

# plot the integrand
fig = plt.figure("Contribution Functions")
ax1 = fig.add_subplot(111)
ax1.plot(τ, C1, label="Linear " * L"S_\nu")
ax1.plot(τ, C2, label="Quadratic " * L"S_\nu")
ax1.set_xlabel(L"\tau_\nu")
ax1.set_ylabel(L"S_\nu(\tau_\nu) \exp(\tau_\nu)")
ax1.set_xlim(0,10)
ax1.legend()
fig.savefig(outdir*"hw3_contribution.pdf")
plt.clf(); plt.close()

# integrate out the contribution functions
nt = 10000           # number of trapezoids
τs = (0.0, 40.0)    # integral bounds for τ
I1 = trap_int(x -> Cν(x, Sν_lin, an...), τs, nt)
I2 = trap_int(x -> Cν(x, Sν_quad, an...), τs, nt)

# compare to E-B approx
println(Sν_lin(1.0, an...))
println(Sν_quad(1.0, an...))

# plot int cont function over source function as a function of a_2
a2 = range(-2.0, 2.0, length=100)
Is = zeros(length(a2))
EB = zeros(length(a2))
for i in 1:length(a2)
    Is[i] = trap_int(x -> Cν(x, Sν_quad, 1.0, 1.0, a2[i]), τs, nt)
    EB[i] = Sν_quad(1.0, 1.0, 1.0, a2[i])
end

# plot the result
fig = plt.figure("Contribution Error")
ax1 = fig.add_subplot(111)
ax1.plot(a2, Is, label="True Intensity")
ax1.plot(a2, EB, ls=":", c="k", label="E-B Intensity")
ax1.set_xlabel("Source Function Quadratic Term " * L"(a_2)")
ax1.set_ylabel(L"I_\nu^+(0, \mu=1)")
ax1.set_xlim(-2.0, 2.0)
ax1.legend()
fig.savefig(outdir*"hw3_quadratic_term.pdf")
plt.clf(); plt.close()

# plot integration error and integral value as function of number of traps
nt = Int.(range(10, 1000, length=100))
τ2 = range(1.0, 100.0, length=100)
It = zeros(length(nt), length(τ2))
for i in 1:length(nt)
    for j in 1:length(τ2)
        It[i,j]= trap_int(x -> Cν(x, Sν_quad, 1.0, 1.0, 1.0), (0.0, τ2[j]), nt[i])
    end
end

# calculate the relative error
re = log10.(abs.(4.0 .- It)./4.0)
extent = [0.0,100,10,1000]

# plot it
fig = plt.figure("Integration Error")#, aspect="equal")
ax1 = fig.add_subplot(111)
img = ax1.imshow(re, cmap="plasma", interpolation="nearest", extent=extent, origin="lower")
cbr = fig.colorbar(img)
cbr.set_label("log(Relative Integration Error)")
ax1.set_xlabel(L"\tau_\nu" * " Upper Bound of Integration")
ax1.set_ylabel("Number of Trapezoids")
ax1.set_aspect("auto")
fig.savefig(outdir*"hw3_integration_error.pdf")
plt.clf(); plt.close()
