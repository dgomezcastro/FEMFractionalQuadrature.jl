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

