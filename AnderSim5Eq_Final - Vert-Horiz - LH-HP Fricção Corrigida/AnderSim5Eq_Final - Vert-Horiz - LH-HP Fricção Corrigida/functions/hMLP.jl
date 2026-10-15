function hMLP(u,Np,nE,V,P)
    # Compute cell averages
    uh = V\u; uh[2:Np,:].=0; uavg = V*uh; v = uavg[1,:];
    
    # max and min cell average per cell in vertex position
    um = zeros(1,nE+1)
    up = zeros(1,nE+1)
    um= vcat(v[1], v); up= vcat(v, v[end]);
    umax = zeros(nE+1)
    umin = zeros(nE+1)
    for i=1:nE+1
        umax[i] = maximum([(um[i]),(up[i])])
        umin[i] = minimum([(um[i]),(up[i])])
    end
    
    # Apply slope limiter as needed
    u_limit = u; eps0=1.0e-8
    
    # find end values of each element
    ue1 = u[1,:]; ue2 = u[end,:];
    
    # find cell averages
    vnE = v; 
    vnEm1 = zeros(nE,1)
    vnEp1 = zeros(nE,1)

    auxvnEm1 = vcat(v[1],v[1:nE-1]); auxvnEp1 = vcat(v[2:nE],v[nE]); 
    vnEm1[:] = auxvnEm1[:]'
    vnEp1[:] = auxvnEp1[:]'
    
    # Apply reconstruction to find elements in need of limiting
    ve1 = zeros(nE,1)
    ve2 = zeros(nE,1)
    for i=1:nE
        ve1[i] = vnE[i] - MINMOD((vnE[i]-ue1[i]),(vnE[i]-vnEm1[i]),(vnEp1[i]-vnE[i]))
        ve2[i] = vnE[i] + MINMOD((ue2[i]-vnE[i]),(vnE[i]-vnEm1[i]),(vnEp1[i]-vnE[i]))
    end
    aux1 = (ve1[:]-ue1[:])./ue1[:]
    aux2 = (ve2[:]-ue2[:])./ue2[:]
    for i=1:nE
        aux1[i] = abs(aux1[i])
        aux2[i] = abs(aux2[i])
    end
    aux3 = zeros(nE,1)
    for i=1:nE
        if (aux1[i]>eps0 || aux2[i]>eps0) == 1.0
            aux3[i] = true
        else
            aux3[i] = false
        end
    end

    ids = zeros(nE)
    contador = 0
    for i=1:nE
        if aux3[i] == 1.0
            contador = contador + 1
            ids[contador] = i 
        end
    end
    ids = Int.(ids[1:contador])
    

    # ChecnE to see if any elements require limiting
    if(!isempty(ids))
        # create piecewise linear solution for limiting on specified elements
        uhl = V\u[:,ids]; uhl[3:Np,:].=0; ul = V*uhl;
        uh0 = V\u[:,ids]; uh0[2:Np,:].=0; u0 = V*uh0;
        # P1j mode
        P1j = ul-u0;
            # apply slope limiter to selected elements
        u_limit[:,ids] = limiter(P1j,v[ids],umax,umin,Np,ids)
    end
    # for i=1:size(u_limit,2)
    #     idx = i;
    #     for j=1:size(u_limit,1) 
    #         if (u_limit[j,i] >= umin[idx] && u_limit[j,i] <= umax[idx])!=1 , 
    #             disp("dmp")
    #         end
    #         idx = idx+1;
    #     end  
    # end
    return u_limit'
end


    function limiter(P1j,s_avg,umax,umin,nN,ids)
    
        ## MLPu limiter
        ver_pos=[1 nN]; mlpu = zeros(1,2); u_limit=zeros(size(P1j)); tol = 1e-8;
        a = P1j
        b = umin
        c = umax
        d = s_avg
        e = ver_pos
        global a, b, c, d, e
        for i=1:size(P1j,2)
            idx = ids[i];
            for p=1:2
                rate=maximum([(umin[idx]-s_avg[i])/(P1j[ver_pos[p],i]+tol),(umax[idx]-s_avg[i])/(P1j[ver_pos[p],i]+tol)]);
            
                if abs(P1j[ver_pos[p],i])>= tol || abs(abs(P1j[ver_pos[p],i])-tol)<eps(Float64)
                    mlpu[p] = minimum([1,rate]); 
                else
                    mlpu[p] = 1;
                end
                idx = idx+1;
            end
            phi_mlpu = minimum(mlpu);
            u_limit[:,i] = s_avg[i].+ phi_mlpu*P1j[:,i];        
        end
        return u_limit
    end

    function minmod2(v)

        # function mfunc = minmod(v)
        # Purpose: Implement the midmod function v is a vector
        
        m = size(v,1); mfunc = zeros(1,size(v,2));
        s = sum(sign(v),1)/m;
        
        ids = find(abs(s)==1);
        if(~isempty(ids))
          mfunc[ids] = s[ids].*minimum([abs(v[:,ids]),[],1]); 
        end
        return mfunc
    end