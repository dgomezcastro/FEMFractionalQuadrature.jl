abstract type AbstractFEM2dBasis <: AbstractFEMBasis end

using Triangulate

include("FEMNeumann.jl")
include("FEMDirichlet.jl")
include("WFEM.jl")
include("kernel_conv.jl")
include("Quad.jl")
include("assemble.jl")

include("particular_geometries/unit_circle.jl")
include("particular_geometries/unit_square.jl")