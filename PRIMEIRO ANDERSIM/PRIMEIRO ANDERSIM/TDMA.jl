function TDMA(B,C,A,D)
    n = length(D)
    x = zeros(n)
    a = A
    b = B
    c = C
    d = D
    # Solver
    c[1] = c[1] / b[1]
    d[1] = d[1] / b[1]
    #
    for i = 2:n
        Aux = b[i] - a[i] * c[i-1]
        c[i] = c[i] / Aux
        d[i] = (d[i] - a[i] * d[i-1]) / Aux
    end
    #
    x[n] = D[n]
    for i = (n - 1):-1:1
        x[i] = d[i] - c[i] * x[i+1]
    end
	return x
end
