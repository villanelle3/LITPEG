function frgrid(LMIN, LMAX, N, P, C)
    # Create mesh
    T = Float64
    E = 1
    vertexPoints = collect(LinRange(LMIN, LMAX, N + 1))
    elementsVertices = Matrix{T}(undef, N, 2)
    for i ∈ 1:length(vertexPoints) - 1
        elementsVertices[i, :] = [vertexPoints[i], vertexPoints[i + 1]]
    end
    mesh = Mesh(E, P, "lobatto", elementsVertices, C; T = T)
    innerPoints = vcat([collect(transpose(mesh.elements[n].innerPoints)) for n ∈ 1:meshSize(mesh)]...)    
    return mesh, innerPoints
end