include("feature_assemble_Stiff.jl")
include("feature_utils.jl")
import .FEM_1d: FEM1d_generate_mesh, FEM1d_connectivity_matrix, FEM1d_linear_elembase
import .Quadrature: generate_quadrature, generate_optimal_quadrature, find_element
import .AssembleStiff: assemble_stiffness
import .MatricesSaveLoad: save_matrix

using LinearAlgebra, YAML, Base.Threads, Plots, Printf, PrettyTables, SpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=1, 
    framestyle=:box, 
    label=nothing, 
    grid=false
)

# Load the configuration file
config_file = "config.yml"
configuration = YAML.load_file(config_file)

# Extract values from the loaded configuration
a = configuration["left_endpoint"]
b = configuration["right_endpoint"]
Ra = configuration["cube_left_endpoint"]
Rb = configuration["cube_right_endpoint"]
h = configuration["mesh_size"]
ρ = configuration["quadrature_size"]
s = configuration["fractional_order"]

α = 2*s
cnst = (2^α * (gamma(1/2 + s))) / (pi^(1/2) * abs(gamma(-s)))

mesh = FEM1d_generate_mesh(a, b, h)
nNode = length(mesh)

quadrature = generate_quadrature(Ra, Rb, ρ)
#quadrature, Ra, Rb = generate_optimal_quadrature(ρ, s)
nQuad = length(quadrature)

println("___FEM and Quadrature on (a,b)=($a ,$b) and QR = ($Ra, $Rb)___ 
mesh size: $h 
quadraure size: $ρ 
fractional order: $s.")


conn = FEM1d_connectivity_matrix(nNode)

@time A = assemble_stiffness(mesh, h, quadrature, ρ, s)

directory = @sprintf("Matrices/s|%.2f_h|%.3f_rho|%.3f_R|%.2f/", s, h, ρ, Rb) #save the unscaled matrix (modify this, it's a residue from when I was uncertain about the constant)
save_matrix(A, directory, s, h, ρ, a, b, Ra, Rb, nNode, nQuad)

A = A * (cnst/2) #scaling the matrix


