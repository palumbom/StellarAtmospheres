# use packages
using Pkg; Pkg.activate(".")
using StellarAtmospheres
using LaTeXStrings
using Plots; pyplot()

# calc blackbody curve in frequency space for a few temperatures
T = [2500.0, 5700.0, 10000.0]
ν = range(1e4, 3e15, length=1000) # in Hz
I = zeros(length(ν), length(T))
for i in 1:length(T)
    I[:,i] = Bν(T[i], ν)
end

# plot it
plot(ν, I[:,1])
plot!(xaxis=:log10,
      yaxis=:log10,
      xlabel="Wavelength (cm)",
      ylabel=L"B_ν",
      legend=:topleft)
