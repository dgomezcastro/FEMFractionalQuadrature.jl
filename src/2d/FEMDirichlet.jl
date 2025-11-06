export PLFEMBasis2dDirichlet

""" 
    Stores a Neumann basis and `interior_point` such that 
    `interior_point[i]` is the index in basisNeumann.mesh.pointlist of the i-th interior point
"""
struct PLFEMBasis2dDirichlet <: AbstractFEM2dBasis
    basisNeumann::PLFEMBasis2dNeumann
    interior_point::Vector{Int64}

    function PLFEMBasis2dDirichlet(mesh::Triangulate.TriangulateIO)
        basisNeumann = PLFEMBasis2dNeumann(mesh)
        is_interior_point = Vector{Bool}(undef, numberofpoints(mesh))
        is_interior_point .= true
        is_interior_point[mesh.segmentlist[1, :]] .= false
        is_interior_point[mesh.segmentlist[2, :]] .= false

        interior_point = zeros(sum(is_interior_point))
        i = 1
        for k in eachindex(interior_point)
            while !is_interior_point[i]
                i = i + 1
            end
            interior_point[k] = i
        end

        return new(basisNeumann, interior_point)
    end
end
dimension(basis::PLFEMBasis2dDirichlet) = length(basis.interior_point)
(basis::PLFEMBasis2dDirichlet)(i, x) = basis.basisNeumann(interior_point[i], x)

mesh(basis::PLFEMBasis2dDirichlet) = basis.basisNeumann.mesh

integral(basis::PLFEMBasis2dDirichlet, i, f::Function) = integral(basis::PLFEMBasis2dNeumann, basis.interior_point[i], f::Function)

@inline function sign_area(p1, p2, p3)
    return (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])
end