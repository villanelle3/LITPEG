function FR(flow, u0, α, ρ, α0, ρ0,vel, mesh, CFL, Δt,Δx,P,θ,dif_art,Ø,LIXO1,LIXO2,fi0,τi,τw,αGF0,D,a1,a2,a3,a4,a5,a6,a7,a8,α1, ρ1, u1,α2, ρ2, u2, wavev, flowpat, Sw, Si) 
    
    #stepRK3SSPS3!(flow, u0, α, ρ, α0, ρ0, Δt,Δx,P,θ,dif_art,vel,Ø, mesh, CFL)
    # stepRK2SSP!(flow, mesh, CFL)
    stepEuler!(flow, u0, α, ρ, α0, ρ0, Δt,Δx,P,θ,dif_art,vel,Ø,fi0,τi,τw,αGF0,D,a1,a2,a3,a4,a5,a6,a7,a8,α1, ρ1, u1,α2, ρ2, u2,wavev, flowpat, Sw, Si, mesh, CFL)
    #aux3 = zeros(size(flow.uδI))
    #aux4 = zeros(Int(size(aux3,1)-3))
    #aux3[:] = flow.uδI[:, 1]
    #aux4 = Centers(flow)
     return (flow.uδD[:,:,1])
end