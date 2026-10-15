function hmlp2(u, Np, nE, V, P)
    # Compute cell averages
    uh = V \ u
    uh[2:Np, :] .= 0
    uavg = V * uh
    v = uavg[1, :]
    
    # max and min cell average per cell in vertex position
    um = [v[1]; v]
    up = [v; v[end]]
    umax = max.(um, up)
    umin = min.(um, up)
    
    # Apply slope limiter as needed
    u_limit = copy(u)
    eps0 = 1.0e-8
    
    # find end values of each element
    ue1 = u[1, :]
    ue2 = u[end, :]
    
    # find cell averages
    vnE = v
    vnEm1 = [v[1]; v[1:nE-1]]
    vnEp1 = [v[2:nE]; v[nE]]
    
    # Apply reconstruction to find elements in need of limiting
    ve1 = vnE' .- minmod2([(vnE .- ue1)'; (vnE .- vnEm1)'; (vnEp1 .- vnE)'])
    ve2 = vnE' .+ minmod2([(ue2 .- vnE)'; (vnE .- vnEm1)'; (vnEp1 .- vnE)'])

    col_mask = any((abs.(ve1 .- ue1) .> eps0) .| (abs.(ve2 .- ue2) .> eps0), dims=1)
    ids = findall(vec(col_mask))

    # Check to see if any elements require limiting
    if !isempty(ids)
        # create piecewise linear solution for limiting on specified elements
        uhl = V \ u[:, ids]
        uhl[3:Np, :] .= 0
        ul = V * uhl
        uh0 = V \ u[:, ids]
        uh0[2:Np, :] .= 0
        u0 = V * uh0
        
        # P1j mode
        P1j = ul - u0
        
        # apply slope limiter to selected elements
        u_limit[:, ids] .= limiter(P1j, v[ids], umax, umin, Np, ids)
    end
    
    return u_limit'
end

function limiter(P1j, s_avg, umax, umin, nN, ids)
    # MLPu limiter
    ver_pos = [1, nN]
    mlpu = zeros(2)
    u_limit = zeros(size(P1j))
    tol = 1e-8
    
    for i in eachindex(ids)
        idx = ids[i]
        for p in 1:2
            rate = max((umin[idx] - s_avg[i]) / (P1j[ver_pos[p], i] + tol),
                       (umax[idx] - s_avg[i]) / (P1j[ver_pos[p], i] + tol))
            
            if abs(P1j[ver_pos[p], i]) >= tol || abs(abs(P1j[ver_pos[p], i]) - tol) < eps()
                mlpu[p] = min(1, rate)
            else
                mlpu[p] = 1
            end
            idx += 1
        end
        
        phi_mlpu = minimum(mlpu)
        u_limit[:, i] .= s_avg[i] .+ phi_mlpu .* P1j[:, i]
    end
    
    return u_limit
end