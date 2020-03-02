"""
Author: Michael Palumbo
Created: April 2019
Contact: mlp95@psu.edu
"""

module StellarAtmospheres

# abbreviate AbstractArray
const AA = AbstractArray
const AF = AbstractFloat

# bring in constants and thermal stuff
include("constants.jl")
include("utils.jl")
include("expint.jl")
include("integrate.jl")
# include("gauss_legendre.jl")
include("differentiate.jl")
include("axis_strings.jl")
include("thermal.jl")
include("spec_tables.jl")
include("partition.jl")
include("opacity.jl")
include("transport.jl")

# export things
export Bν, Bλ, Bν̃, λmax, νmax,
       λ2ν, ν2λ, λ2ν̃, ν̃2λ, ν̃2ν,
       ν2ν̃, SνPoly, SνPlanck,
       Cν, Iν₀, IνEB, ℱν₀,
       ℱντ, ℱτ, ℱνEB, Hν₀,
       HνEB, SB, deriv, deriv2

end # module
