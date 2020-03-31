using Pkg; Pkg.activate(".")
using Revise
using Dierckx
using BenchmarkTools
using StellarAtmospheres; SA = StellarAtmospheres;
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib;

# set up plot output
outdir = "/Users/michael/Desktop/ASTRO530/figures/"
mpl.style.use("atmospheres.mplstyle"); plt.ioff()

# get VALIIIC atmosphere
dfv = SA.dfv

# get stuff at τ=1 surface
ind = 49
temp = dfv.T[ind]
Pg = dfv.Pg_Ptot[ind] .* dfv.Ptot[ind]
Pe1 = SA.P_from_nkt.(dfv.n_e[ind], temp)
Pe = SA.calc_Pe.(temp, Pg, Pe1)
ξ = dfv.V[49] * 1000.0 * 100 # convert km/s -> cm/s
N_NE = 0.8609 # ground fraction
nH = dfv.n_H[ind]
ρ = dfv.ρ[ind]
f_ground = 0.8609
f_neutral = 6.5418e-4
m = 3.817541e-23

# make LineParams object instances for NaD lines
NaD2 = LineParams(element="Na", n=3, λ₀=5890.0, A=[1e8*6.16e-1/(4π)], m=m, gu=4, gl=2, logC4=-15.17)
NaD1 = LineParams(element="Na", n=3, λ₀=5896.0, A=[1e8*6.14e-1/(4π)], m=m, gu=2, gl=2, logC4=-15.33)

# compare some values
SE = SA.calc_SE(5890.0, temp)
# println("γn = " * string(SA.calc_γn(NaD2)))
# println("log10(γ4) = " * string(log10(SA.calc_γ4(Pe, temp, NaD2))))
# println("log10(γ6) = " * string(log10(SA.calc_γ6(Pg, temp, NaD2))))
# println("ΔλD = " * string(SA.calc_ΔλD(temp, ξ, NaD2)))
# println("SE factor = " * string(SE))

# get alpha and convert to opacity
alph = SA.calc_α(5890.0, Pe, Pg, temp, ξ, NaD2)
opac = alph * f_ground * f_neutral * SA.abundance_for_element("Na") * nH / ρ

# waves = range(5880.0, 5900.0, length=1000)
# alpha = SA.calc_α(waves, Pe, Pg, temp, ξ, NaD2)
# plt.plot(waves, alpha)

# # continuous opacity
# PT = (1.0 + SA.ΦT(temp, "H")/Pe)
# pg = SA.sum_abundance_weights()
# κcont = SA.κ_tot(5890.0, temp, Pe, Pg)
# κHmbf = SA.κ_H_minus_bf(5890.0, temp, Pe) * SE / PT / pg
# κHmff = SA.κ_H_minus_ff(5890.0, temp, Pe) / PT / pg
# κHbf = SA.κ_H_bf(5890.0, temp, Pe) * SE / PT / pg
# κHff = SA.κ_H_ff(5890.0, temp, Pe) * SE / PT / pg
# κe = SA.κ_e(Pe, Pg) / pg

# # print them
# @show κcont
# @show κHmbf
# @show κHmff
# @show κHbf
# @show κHff
# @show κe
