
module MatricesSaveLoad

using Printf
using DelimitedFiles
using FilePathsBase

function save_matrix(A, directory, s, h, ρ, a, b, Ra, Rb, nNode, nQuad)
    # Create directory if it doesn't exist
    if !isdir(directory)
        mkdir(directory)
    end

    # Write input data to a file
    open(joinpath(directory, "input_data.txt"), "w") do file
        write(file, @sprintf("Dirichlet problem for Fractional Laplacian with fractional order s=%.2f on (a,b)=(%.2f, %.2f) \n quadrature for (Ra, Rb)= (%.2f, %.2f) \n Stiffness matrix for : \n rho: %.3f \n h: %.3f \n number of nodes: %d \n number of quadrature points: %d",
                            s, a, b, Ra, Rb, ρ, h, nNode, nQuad))
    end

    # Save matrix A to file
    filename = @sprintf("Stiff_s|%.2f_h|%.3f_rho|%.3f.txt", s, h, ρ)
    writedlm(joinpath(directory, filename), A)
end

function load_matrix(directory)
    return readdlm(directory)
end

end

module LpNorms

function L1_trapz(func, h)

    return sum(((func[i] + func[i+1]) * h/ 2) for i in range(1, length(func)-1))

end

function L2_trapz(func, h)

    return sqrt(sum(((func[i]^2 + func[i+1]^2) * h / 2) for i in range(1, length(func)-1)))

end

function L2norm(rho2, rho, Dx)
    norm = 0.0
    N = length(rho)
    
    norm += ((rho[1] - rho2[1])^2 + 2 * ((rho[1] + rho[2]) / 2 - rho2[2])^2 +
             (rho[2]-rho2[3])^2) * Dx / 2

    for i in 2:(N-2)
        norm += ((rho[i] - rho2[2*i-1])^2 + 2 * ((rho[i] + rho[i+1]) / 2 - rho2[2 * i])^2 +
             (rho[i+1]-rho2[2*i+1])^2) * Dx / 2
    end

    norm += ((rho[N-1] - rho2[2*(N-1)-1])^2 + 2 * ((rho[N-1] + rho[N]) / 2 - rho2[2*(N -1)])^2 +
             (rho[N]-rho2[2*(N-1)+1])^2) * Dx / 2

    return sqrt(norm)
end

function L2norm_Simpson(rho2, rho, Dx)
    norm = 0.0
    N = length(rho)
    
    norm += ((rho[1] - rho2[1])^2 + 2 * 2*((rho[1] + rho[2]) / 2 - rho2[2])^2 +
            4 * (((rho[1] + (rho[1] + rho[2]) / 2) / 2) - (rho2[1]+ rho2[2])/2)^2 +
            4 * (((rho[2] + (rho[1] + rho[2]) / 2) / 2) - (rho2[2]+ rho2[3])/2)^2 + 
             (rho[2]-rho2[3])^2) * Dx / 6

    for i in 2:(N-2)
        norm += ((rho[i] - rho2[2*i-1])^2 + 2 * ((rho[i] + rho[i+1]) / 2 - rho2[2 * i])^2 +
                4 * (((rho[i] + (rho[i] + rho[i+1]) / 2) / 2) - (rho2[2*i-1]+ rho2[2 * i])/2)^2 +
                4 * (((rho[i+1] + (rho[i] + rho[i+1]) / 2) / 2) - (rho2[2*i]+ rho2[2*i+1])/2)^2 + 
             (rho[i+1]-rho2[2*i+1])^2) * Dx / 6
    end

    norm += ((rho[N-1] - rho2[2*(N-1)-1])^2 + 2 * ((rho[N-1] + rho[N]) / 2 - rho2[2*(N -1)])^2 +
            4 * (((rho[N-1] + (rho[N-1] + rho[N]) / 2) / 2) - (rho2[2*(N-1)-1]+ rho2[2*(N-1)])/2)^2 +
            4 * (((rho[N] + (rho[N-1] + rho[N]) / 2) / 2) - (rho2[2*(N-1)]+ rho2[2*(N-1)+1])/2)^2 + 
             (rho[N]-rho2[2*(N-1)+1])^2) * Dx / 6

    return sqrt(norm)
end

end
