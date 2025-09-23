f(x) = 1.0
F = h*[f(x) * (dist_boundary_regular_UnitInterval(x)^s) for x in mesh[2:nNode-1]]
F = [0; F; 0]

