using LinearAlgebra, Base.Threads
function fractional_Laplacian_stiffness_matrix(
    s::Float64,
    fe::FEM1d,
    quad::QuadratureParameters,
)
    N = nbasisfunctions(fe)
    A = zeros(N, N)

    @threads for i = 1:N
        for j = i:N
            A[i, j] = sobolevslobodeckiiprod_basis(s, fe, quad, i, j)
        end
    end
    return Symmetric(A)
end

function rhs(f::Function, fe::FEM1d; Δx = 1e-5)
    x = fe.mesh[1]:Δx:fe.mesh[end]
    return [Δx * sum(f.(x) .* basis(fe, i).(x)) for i = 1:nbasisfunctions(fe)]
end
