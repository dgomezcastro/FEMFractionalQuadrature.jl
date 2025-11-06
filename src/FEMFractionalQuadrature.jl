module FEMFractionalQuadrature

using LinearAlgebra, Base.Threads, Logging

include("EpsteinLib.jl") # TODO: Replace by `import EpsteinLib` when package is available. Update Project.toml

include("1d/1d.jl")
include("2d/2d.jl")
include("ExtendableFEM/ExtendableWFEM.jl")

end # module FEMFractionalQuadrature##