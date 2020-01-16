"""
Author: Michael Palumbo
Created: April 2019
Contact: mlp95@psu.edu
"""

module StellarAtmospheres

# bring in constants and thermal stuff
include("constants.jl")
include("thermal.jl")
include("pretty_plots.jl")

# export things
export Bν, Bλ, λmax, νmax, λ2ν, ν2λ

end # module
