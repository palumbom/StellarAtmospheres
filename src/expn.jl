using Conda; ENV["CONDA_JL_HOME"] = "/usr/local/anaconda3/envs/conda_jl";
using PyCall; py_expn = pyimport("scipy.special").expn;

function expn(n::Int, x::T) where T<:Real
    return py_expn(n, x)
end
