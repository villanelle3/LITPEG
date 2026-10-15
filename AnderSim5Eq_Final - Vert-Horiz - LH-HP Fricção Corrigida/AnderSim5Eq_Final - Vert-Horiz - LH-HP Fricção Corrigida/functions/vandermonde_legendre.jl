function vandermonde_legendre(N::Int, r::Vector{Float64})
    # Pré-aloca a matriz
    V = zeros(length(r), N+1)
    
    # Preenche a primeira coluna com P0(x) = 1
    V[:, 1] .= 1.0
    
    # Caso N ≥ 1, preenche a segunda coluna com P1(x) = x
    if N ≥ 1
        V[:, 2] = r
    end
    
    # Usa a relação de recursão para n ≥ 2
    for n in 1:N-1
        V[:, n+2] = ((2n + 1) .* r .* V[:, n+1] - n * V[:, n]) / (n + 1)
    end
    
    return V
end