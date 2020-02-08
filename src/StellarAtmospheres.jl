"""
Author: Michael Palumbo
Created: April 2019
Contact: mlp95@psu.edu
"""

module StellarAtmospheres

# bring in constants and thermal stuff
include("constants.jl")
include("utils.jl")
include("expint.jl")
include("integrate.jl")
include("axis_strings.jl")
include("thermal.jl")
include("transport.jl")

# export things
export Bν, Bλ, Bν̃, λmax, νmax,
       λ2ν, ν2λ, λ2ν̃, ν̃2λ, ν̃2ν,
       ν2ν̃, trap_int, Sν, Cν,
       ℱν₀, Hν₀

end # module
