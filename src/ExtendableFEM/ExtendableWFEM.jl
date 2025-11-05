# TODO: Create tests

export ExtendableWFEMBasis

struct ExtendableWFEMBasis
    s::Real
    distance::Function
    quotientspace::ExtendableFEM.FESpace
    PEs::Vector{ExtendableFEM.PointEvaluator}

    function WFEMBasis(s, quotientspace::ExtendableFEM.FESpace, distance::Function)

        φs = [FEVector(quotientspace) for i in 1:quotientspace.ndofs]
        for i = 1:quotientspace.ndofs
            φs[i].entries[i] = 1.0
        end
        PEs = [PointEvaluator([(1, Identity)], FEVector(φ)) for φ in φs]

        return new(s, quotienbasis, PEs, distance)
    end
end

dimension(basis::ExtendableWFEMBasis) = size(basis.Nodes, 2)

"""
function evaluating the ϕ_i basis function at the point P
"""
function ϕ(basis::ExtendableWFEMBasis, i, P)
    result = zeros(2)
    evaluate!(result, basis.PEs[i], P)
    return basis.distance(P)^s * result[1]
end

# Adapted from https://wias-pdelib.github.io/ExtendableFEMBase.jl/dev/module_examples/Example200_LowLevelPoisson/
function integral(basis::ExtendableWFEMBasis, f::Function)
    fe_space = basis.quotientspace
    xgrid = fe_space.xgrid
    EG = xgrid[UniqueCellGeometries][1]
    FEType = eltype(fe_space)
    L2G = L2GTransformer(EG, xgrid, ON_CELLS)

    qf = QuadratureRule{Float64,EG}(2 * (get_polynomialorder(FEType, EG) - 1))
    weights::Vector{Float64} = qf.w
    xref::Vector{Vector{Float64}} = qf.xref
    nweights::Int = length(weights)
    cellvolumes = xgrid[CellVolumes]

    FEBasis_id = FEEvaluator(fe_space, Identity, qf)
    idvals = FEBasis_id.cvals

    x::Vector{Float64} = zeros(Float64, 2)

    for cell in 1:ncells
        for j in 1:ndofs4cell
            # right-hand side
            temp = 0
            for qp in 1:nweights
                # get global x for quadrature point
                eval_trafo!(x, L2G, xref[qp])
                # evaluate (f(x), v_j(x))
                temp += weights[qp] * idvals[1, j, qp] * f(x)
            end
            # write into global vector
            integral += temp * cellvolumes[cell]
        end
    end

    return integral
end

function integral(basis::WFEMBasis, f::Function)

end