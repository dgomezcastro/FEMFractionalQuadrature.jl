module FEMFractionalQuadrature

abstract type FEM1d end
include("lfem.jl")
include("wfem.jl")
include("quadrature.jl")
include("assemble.jl")

include("exactsolution.jl")

end # module FEMFractionalQuadrature
