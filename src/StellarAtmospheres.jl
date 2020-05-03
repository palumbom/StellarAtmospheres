"""
Author: Michael Palumbo
Created: April 2019
Contact: mlp95@psu.edu
"""

module StellarAtmospheres

# abbreviate AbstractArray
const AA = AbstractArray
const AF = AbstractFloat

# some modules for global use
using RecursiveArrayTools

# bring in constants and thermal stuff
include("constants.jl")
include("LineParams.jl")
include("utils.jl")
include("expint.jl")
include("voigt.jl")
include("integrate.jl")
# include("gauss_legendre.jl")
include("differentiate.jl")
include("axis_strings.jl")
include("partition.jl")
include("thermal.jl")
include("data_tables.jl")
include("abundances.jl")
include("ionizations.jl")
include("nlte.jl")
include("opacity.jl")
include("pressure.jl")
include("transport.jl")

# export things
export Bν, Bλ, Bν̃, λmax, νmax,
       λ2ν, ν2λ, λ2ν̃, ν̃2λ, ν̃2ν,
       ν2ν̃, SνPoly, SνPlanck,
       Cν, Iν₀, IνEB, ℱν₀,
       ℱντ, ℱτ, ℱνEB, Hν₀,
       HνEB, SB, LineParams,
       κ_tot, κ_line, asum,
       elav, ℱν₀_line, expint,
       trap_int, calc_τν

end # module
