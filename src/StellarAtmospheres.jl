"""
Author: Michael Palumbo
Created: April 2019
Contact: mlp95@psu.edu
"""

module StellarAtmospheres

# bring in constants
include("constants.jl")

# define planck functions
Bν(ν::t, T::t) where t<:Real = 2.0*h*ν^3/c^2 * 1.0/(exp(h*ν/(kB*T)) - 1.0)
Bλ(λ::t, T::t) where t<:Real = 2.0*h*c^2/λ^5 * 1.0/(exp(h*c/(λ*kB*T)) - 1.0)

end # module
