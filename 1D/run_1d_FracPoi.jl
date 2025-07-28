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
cnst_2 = (α * 2^(α-1) * gamma(1/2 + s)) / (pi^(1/2) * gamma(1-s))
cnst_3 = s*2^(-1+α)*gamma(1+s)/(pi^2*gamma(1-s))
cnst_4 = (s * 2^(α) * gamma(1/2 + s)) / (pi^(1/2) * gamma(1-s))
cnst_book1 = 1/pi * gamma(1+α) * sin((α*pi)/2)
cnst_book2 = α * (2^(α-1)/sqrt(pi)) * gamma(1/2 + s) / gamma(1-s)
@show cnst, cnst_2, cnst_3, cnst_4
cnst_sol = (pi^(1/2) * 2^(-α)) / (gamma(1/2 + s) * gamma(1+s))
cnst_sol_2 = (2^(1-α) * gamma(1/2)) / (α * gamma(1/2 + s) * gamma(s))
cnst_sol_3 = (gamma(1/2)) / (2^α * gamma(1/2 + s) * gamma(1+s))
sol(x) = cnst_sol * (-x.^2 .+ 1).^s
@show cnst_sol, cnst_sol_2, cnst_sol_3
L1norm_sol = sum(cnst_sol * (1 .- LinRange(a, b, 4097).^2).^s) * 2 / 4096
Mass = (pi) / (2^α * gamma(s+1/2) * gamma(s+3/2))
@show L1norm_sol, Mass

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

@time A = assemble_stiffness(mesh, h, quadrature, ρ, s, 1)

directory = @sprintf("Matrices/s|%.2f_h|%.3f_rho|%.3f_R|%.2f/", s, h, ρ, Rb)
save_matrix(A, directory, s, h, ρ, a, b, Ra, Rb, nNode, nQuad)

A = A * (cnst/2)

f(x) = 1.0
F = h*[f(x) for x in mesh[2:nNode-1]]
F = [0; F; 0]

U = (A \ F)

U_scal = L1norm_sol/(sum(U)*h) * U

@show sum(U)*h
@show L1norm_sol / cnst_sol
@show L1norm_sol/(sum(U)*h)
@show sol(mesh)[Int(nNode ÷ 2)] / U[Int(nNode ÷ 2)]
@show 1/(2^s)

#formatted_figname = @sprintf("Figures/FEFrac_s|%.2f_h|%.3f_rho|%.3f.png", s, h, ρ)
#cnst = U[Int(nNode ÷ 2)] / sol(mesh)[Int(nNode ÷ 2)]
plot(mesh, U; label="\$u_h\$", marker=:circle, show=true, size = (800,600),
tickfontsize = 10, legendfontsize=15, dpi=100)
plot!(mesh, sol(mesh), label="\$u^{*}\$", lw=1.5, dpi=200)
#savefig(formatted_figname) 
readline()


