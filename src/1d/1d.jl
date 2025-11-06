abstract type AbstractFEM1dBasis <: AbstractFEMBasis end

include("FEM.jl")
include("WFEM.jl")
include("Quad.jl")
include("assemble.jl")
include("problem.jl")
include("L2norm.jl")