


function generate_mesh_UnitSquare(h::Float64)

    hh = h^2 / 2 #area of each triangle

    triin = Triangulate.TriangulateIO()
    triin.pointlist = Matrix{Cdouble}([0.0 0.0; 1.0 0.0; 1.0 0.5; 0.0 0.5; 0.0 1.0; 1.0 1.0]')
    triin.segmentlist = Matrix{Cint}([1 2; 2 3; 3 4; 4 5; 5 6; 6 3; 4 1]')
    triin.segmentmarkerlist = Vector{Int32}([1, 2, 3, 4, 5, 6, 7])

    maxarea = hh
    area = @sprintf("%.15f", maxarea)

    (triout, vorout) = triangulate("pa$(area)DQ", triin)
    return triout

end


function WFEM2d_generate_mesh_Square(h::Float64, b::Float64)

    n = Int(2b / h)

    xs = collect(range(-b, b, length=n + 1))
    ys = collect(range(-b, b, length=n + 1))

    points = [(x, y) for y in ys for x in xs]

    node(i, j) = j * (n + 1) + i + 1

    tris = Vector{NTuple{3,Int}}()
    for j in 0:n-1, i in 0:n-1
        p1 = node(i, j)
        p2 = node(i + 1, j)
        p3 = node(i + 1, j + 1)
        p4 = node(i, j + 1)

        push!(tris, (p1, p2, p3))
        push!(tris, (p1, p3, p4))
    end

    triout = Triangulate.TriangulateIO()
    triout.pointlist = reshape(vcat([Float64[x, y] for (x, y) in points]...), 2, :)
    triout.trianglelist = reshape(vcat([Int32[a, b, c] for (a, b, c) in tris]...), 3, :)

    return triout
end

function WFEMBasis2dDirichletUnitCircle_Square(h::Float64, s::Float64; δ::Function=P -> max(1 - norm(P)^2, 0.0))
    mesh = WFEM2d_generate_mesh_Square(h, 1.0)
    return WFEMBasis2dDirichlet(s, mesh, δ)
end

function nodes_in_unit_circle(triout::Triangulate.TriangulateIO)
    nodes = triout.pointlist
    return vec(sum(nodes .^ 2, dims=1) .< 1.0)
end

function nodes_triangles_intersection_unit_circle(triout::Triangulate.TriangulateIO)
    nodes = triout.pointlist
    triangles = triout.trianglelist

    n_nodes = size(nodes, 2)
    n_triangles = size(triangles, 2)

    inside = vec(sum(nodes .^ 2, dims=1) .< 1.0)

    node_mask = falses(n_nodes)

    for t in 1:n_triangles
        i, j, k = triangles[:, t]

        if inside[i] || inside[j] || inside[k]
            node_mask[i] = true
            node_mask[j] = true
            node_mask[k] = true
        end
    end

    return node_mask
end



