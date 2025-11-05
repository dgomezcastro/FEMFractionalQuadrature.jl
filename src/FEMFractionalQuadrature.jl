module FEMFractionalQuadrature

using LinearAlgebra, Base.Threads, Logging

include("1d/1d.jl")
include("2d/2d.jl")
include("ExtendableFEM/ExtendableWFEM.jl")

end # module FEMFractionalQuadrature##