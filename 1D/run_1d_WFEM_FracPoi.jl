include("feature_WFEM_1d.jl")
include("feature_WFEM_assemble_Stiff.jl")
include("feature_utils.jl")
import .WFEM_1d: WFEM1d_generate_mesh, WFEM1d_connectivity_matrix, WFEM1d_linear_elembase, dist_boundary, WFEM1d_ϕ
import .Quadrature: generate_quadrature, generate_optimal_quadrature, find_element
import .AssembleStiff_WFEM: assemble_stiffness_WFEM
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

mesh = WFEM1d_generate_mesh(a, b, h)
nNode = length(mesh)

Dist_Boundary = [dist_boundary(mesh[i], a, b) for i in 1:nNode]

plot(mesh, Dist_Boundary.^s, show=true)
sleep(0.01)

quadrature = generate_quadrature(Ra, Rb, ρ)
#quadrature, Ra, Rb = generate_optimal_quadrature(ρ, s)
nQuad = length(quadrature)

println("___FEM and Quadrature on (a,b)=($a ,$b) and QR = ($Ra, $Rb)___ 
mesh size: $h 
quadraure size: $ρ 
fractional order: $s.")


conn = WFEM1d_connectivity_matrix(nNode)

@time A = assemble_stiffness_WFEM(mesh, h, a, b, quadrature, ρ, s, 1)

directory = @sprintf("Matrices_WFEM/s|%.2f_h|%.3f_rho|%.4f_R|%.2f/", s, h, ρ, Rb)
save_matrix(A, directory, s, h, ρ, a, b, Ra, Rb, nNode, nQuad)

A = A * (cnst/2)

f(x) = 1.0
F = h*[f(x) * (dist_boundary(x, a, b)^s) for x in mesh[2:nNode-1]]
F = [0; F; 0]

U = (A \ F)

U_scal = L1norm_sol/(sum(U)*h) * U

@show sum(U)*h
@show L1norm_sol / cnst_sol
@show L1norm_sol/(sum(U)*h)
@show sol(mesh[2:nNode-1])[Int(nNode ÷ 2)] / U[Int(nNode ÷ 2)]
@show 1/(2^s)

xx = collect(a-h/4:h/4:b+h/4) #finer mesh
NN = length(xx)
global U_WFEM = zeros(length(xx))
for i in 2:nNode-1
    global U_WFEM = U_WFEM .+ U[i] * [WFEM1d_ϕ(mesh[i], xx[j], h, a, b, s) for j in 1:NN]
end

#formatted_figname = @sprintf("Figures/FEFrac_s|%.2f_h|%.3f_rho|%.3f.png", s, h, ρ)
#cnst = U[Int(nNode ÷ 2)] / sol(mesh)[Int(nNode ÷ 2)]
plot(mesh, U; label="\$(u/d^s)_h\$", marker=:circle, show=true, size = (800,600),
tickfontsize = 10, legendfontsize=15, dpi=100)
plot!(mesh[2:nNode-1], sol(mesh[2:nNode-1]) ./ (Dist_Boundary[2:nNode-1].^s), label="\$(u^{*}/d^s)\$", lw=1.5, dpi=200)
#savefig(formatted_figname) 
readline()


plot(mesh, U .* (Dist_Boundary.^s); label="\$u_h\$", marker=:circle, show=true, size = (800,600),
tickfontsize = 10, legendfontsize=15, dpi=100)
plot!(mesh[2:nNode-1], sol(mesh[2:nNode-1]), label="\$u^{*}\$", lw=1.5, dpi=200)
#savefig(formatted_figname) 
readline()

plot(xx, U_WFEM; label="\$u_h\$", marker=:circle, show=true, size = (800,600),
tickfontsize = 10, legendfontsize=15, dpi=100)
plot!(xx[2:NN-1], sol(xx[2:NN-1]), label="\$u^{*}\$", lw=1.5, dpi=200)
plot!(mesh, U .* (Dist_Boundary.^s); label="\$u_h\$", marker=:circle, show=true, size = (800,600),
tickfontsize = 10, legendfontsize=15, dpi=100)
#savefig(formatted_figname) 
readline()

distfine = norm(U_WFEM[2:NN-1] - sol(xx[2:NN-1])) * sqrt(h/4)
distst = norm((U .* (Dist_Boundary.^s))[2:nNode-1] - sol(xx[2:nNode-1])) * sqrt(h)
print("Distance on finer mesh: $distfine 
Distance on normal mesh: $distst")
