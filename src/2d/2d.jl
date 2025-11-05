abstract type AbstractFEM2dBasis end

using Triangulate

include("WFEM.jl")
include("kernel_conv.jl")
include("Quad.jl")
include("assemble.jl")

include("particular_geometries/unit_circle.jl")
include("particular_geometries/unit_square.jl")