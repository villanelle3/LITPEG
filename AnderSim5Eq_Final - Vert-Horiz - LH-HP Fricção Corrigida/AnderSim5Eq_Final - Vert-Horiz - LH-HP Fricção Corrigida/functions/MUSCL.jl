function MUSCL(flow, u0, αF, ρF, αF0, ρF0, u, mesh, CFL, Δt, Δx, p, θ, dif_art, uc, α, ρ, fi0, τi, τw, αGF0, D, flowpat, Sw, Si)
  n = length(u)
  U = zeros(n, 1) # pra guardar u no tempo n+1
  g = 9.81
  if flow.PDE == 1
    for i = 2:(n-1)
      gradp = (p[i] - p[i-1]) / Δx
      #=
      F_R = α[i]*ρ[i]*(uL[i]*(sign((uL[i]+uL[i+1])/2)+1)/2 + uL[i+1]*(-sign((uL[i]+uL[i+1])/2)+1)/2)^2
      F_L = α[i-1]*ρ[i-1]*(uL[i-1]*(sign((uL[i-1]+uL[i])/2)+1)/2 + uL[i]*(-sign((uL[i-1]+uL[i])/2)+1)/2)^2
      U[i] = u0[i] - (Δt/Δx)*(F_R-F_L) - Δt*(αF[i]*gradp + αF[i]*ρF[i]*sin(θ) - dif_art[i])
      =#
      if i == 2
        r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
        r_i_menos = 0
        r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
        uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
        uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
        uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
      elseif i == (n - 1)
        r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
        r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
        r_i_mais = 0
        uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        uRmais = u[i+1]
        uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
        uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
      else
        r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
        r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
        r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
        uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
        uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
        uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
      end
      #
      #= Fluxo na face da direita (i+1/2) - (LLF)
        a1 = abs((uLmais^2-u[i]^2)/(uLmais-u[i]+eps()))
        a2 = abs((uRmais^2-uLmais^2)/(uRmais-uLmais+eps()))
        a3 = abs((u[i+1]^2-uRmais^2)/(u[i+1]-uRmais+eps()))
        a = 1 #maximum([a1,a2,a3]) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      =#

      F_mais = uLmais * (sign((u[i] + u[i+1]) / 2) + 1) / 2 + uRmais * (-sign((u[i] + u[i+1]) / 2) + 1) / 2 # ROE
      #F_mais = (uLmais+uRmais)/2 - a*(uRmais-uLmais)/2 # LLF

      # Fluxo na face da esquerda(i-1/2) - (LLF)
      a1 = abs((uLmenos^2 - u[i-1]^2) / (uLmenos - u[i-1] + eps()))
      a2 = abs((uRmenos^2 - uLmenos^2) / (uRmenos - uLmenos + eps()))
      a3 = abs((u[i]^2 - uRmenos^2) / (u[i] - uRmenos + eps()))
      a = 1 #maximum([a1,a2,a3]) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #

      F_menos = uLmenos * (sign((u[i-1] + u[i]) / 2) + 1) / 2 + uRmenos * (-sign((u[i-1] + u[i]) / 2) + 1) / 2# ROE
      #F_menos = (uLmenos+uRmenos)/2 - a*(uRmenos-uLmenos)/2 # LLF
      if flowpat == 1 # ANULAR
        if flow.BC == 1
          U[i] = fi0[i] - (Δt / Δx) * u[i] * (F_mais - F_menos) - Δt * (αF0[i] * gradp - αF0[i] * ρF0[i] * g * sin(θ) + dif_art[i]) / (αF0[i] * ρF0[i]) + Δt * ( (-4) * τw[i-1] / D + 2 * τi[i-1] * (αGF0[i])^(1 / 2) / (D / 2)) / (αF0[i] * ρF0[i])
        else
          U[i] = fi0[i] - (Δt / Δx) * u[i] * (F_mais - F_menos) - Δt * (αF0[i] * gradp - αF0[i] * ρF0[i] * g * sin(θ) + dif_art[i]) / (αF0[i] * ρF0[i]) - Δt * (2 * αGF0[i] * τi[i-1] / ((αGF0[i])^(1 / 2) * D / 2)) / (αF0[i] * ρF0[i])
        end
      elseif flowpat == 2 # ESTRATIFICADO
        if flow.BC == 1
          U[i] = fi0[i] - (Δt / Δx) * u[i] * (F_mais - F_menos) - Δt * (αF0[i] * gradp - αF0[i] * ρF0[i] * g * sin(θ) + dif_art[i]) / (αF0[i] * ρF0[i]) + Δt * ( (-1)*τw[i-1]*Sw[i-1] / (3.14*D^2/4) + τi[i-1] * Si[i-1] / (3.14*D^2/4)) / (αF0[i] * ρF0[i])
        else
          U[i] = fi0[i] - (Δt / Δx) * u[i] * (F_mais - F_menos) - Δt * (αF0[i] * gradp - αF0[i] * ρF0[i] * g * sin(θ) + dif_art[i]) / (αF0[i] * ρF0[i]) + Δt * ( (-1)*τw[i-1]*Sw[i-1] / (3.14*D^2/4) - τi[i-1] * Si[i-1] / (3.14*D^2/4)) / (αF0[i] * ρF0[i])
        end
      end
    end
    #
  elseif flow.PDE == 2
    n -= 1
    vel = u
    u = α
    for i = 2:(n-1)
      if i == 2
        r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
        r_i_menos = 0
        r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
        uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
        uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
        uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
      elseif i == (n - 1)
        r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
        r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
        r_i_mais = 0
        uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        uRmais = u[i+1]
        uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
        uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
      else
        r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
        r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
        r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
        uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
        uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
        uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
        uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
      end
      # Fluxo na face da direita (i+1/2) - (LLF)
      a1 = abs((uLmais - u[i]) / (uLmais - u[i] + eps()))
      a2 = abs((uRmais - uLmais) / (uRmais - uLmais + eps()))
      a3 = abs((u[i+1] - uRmais) / (u[i+1] - uRmais + eps()))
      a = 1 #maximum([a1,a2,a3])
      F_mais = vel[i+1] * (uLmais * (sign(vel[i+1]) + 1) / 2 + uRmais * (-sign(vel[i+1]) + 1) / 2) # ROE
      #F_mais = vel[i+1]*(uLmais+uRmais)/2 - a*(uRmais-uLmais)/2 # LLF

      # Fluxo na face da esquerda(i-1/2) - (LLF)
      a1 = abs((uLmenos - u[i-1]) / (uLmenos - u[i-1] + eps()))
      a2 = abs((uRmenos - uLmenos) / (uRmenos - uLmenos + eps()))
      a3 = abs((u[i] - uRmenos) / (u[i] - uRmenos + eps()))
      a = 1 #maximum([a1,a2,a3])
      F_menos = vel[i] * (uLmenos * (sign(vel[i]) + 1) / 2 + uRmenos * (-sign(vel[i]) + 1) / 2) # ROE
      #F_menos = vel[i]*(uLmenos+uRmenos)/2 - a*(uRmenos-uLmenos)/2 # LLF

      U[i] = fi0[i] - (Δt / Δx) * (F_mais - F_menos)
    end
    #

    #= implementar aqui a equação de energia 
    elseif flow.PDE == 3 
      n -= 1    
    for i in 1:n
        num = αL0[i] * ρL0[i] * uL0[i] + αG0[i] * ρG0[i] * uG0[i]
        den = αL0[i] * ρL0[i] + αG0[i] * ρG0[i]
        vel[i] = num / den
    end
      # Variável conservada: energia total da mistura (com h2 para ambas as fases)
          # Termos fontes
        for i = 2:N
        c4[i] = (P[i] - PØ[i])/Δt                  # -> foi
        c5[i] = g * sin(θ) * (αL0[i] * ρL0[i] * uL0[i] + αG0[i] * ρG0[i] * uG0[i]) # -> foi
        cn = hL0[i]/Cp_L
    if  cn >= T_ext
        c6[i] = ( hL0[i] * S) / (Cp_L * R_tot[i] * A)     # -> foi

        c7[i] = (-T_ext  * S) / (R_tot[i] * A)           #-> foi
    else
        c6[i] = (-hL0[i] * S) / (Cp_L * R_tot[i] * A)    #-> foi 

        c7[i] = (T_ext * S) / (R_tot[i] * A)            #-> foi

    end                                                  
      for i in 1:n
        u[i] = αG0[i] * ρG0[i] * (Λ * hL0[i] + 0.5 * uG0[i]^2) + αL0[i] * ρL0[i] * (hL0[i] + 0.5 * uL0[i]^2)
    end

     # Fluxos MUSCL + Upwind

             for i=2:(n-1)
              if i==2
                   r_i = (u[i]-u[i-1])/(u[i+1]-u[i]+ eps())
                   r_i_menos = 0
                   r_i_mais = (u[i+1]-u[i])/(u[i+2]-u[i+1]+ eps())
                   uLmais = u[i] + 0.5*Limiter(r_i)*(u[i+1]-u[i])
                   uRmais = u[i+1] - 0.5*Limiter(r_i_mais)*(u[i+2]-u[i+1])
                   uLmenos = u[i-1] + 0.5*Limiter(r_i_menos)*(u[i]-u[i-1])
                   uRmenos = u[i] - 0.5*Limiter(r_i)*(u[i+1]-u[i])
                 elseif i==(n-1)
                   r_i = (u[i]-u[i-1])/(u[i+1]-u[i]+ eps())
                   r_i_menos = (u[i-1]-u[i-2])/(u[i]-u[i-1]+eps())
                   r_i_mais = 0
                   uLmais = u[i] + 0.5*Limiter(r_i)*(u[i+1]-u[i])
                   uRmais = u[i+1] 
                   uLmenos = u[i-1] + 0.5*Limiter(r_i_menos)*(u[i]-u[i-1])
                   uRmenos = u[i] - 0.5*Limiter(r_i)*(u[i+1]-u[i])
                 else
                   r_i = (u[i]-u[i-1])/(u[i+1]-u[i]+eps())
                   r_i_menos = (u[i-1]-u[i-2])/(u[i]-u[i-1]+eps())
                   r_i_mais = (u[i+1]-u[i])/(u[i+2]-u[i+1]+eps())
                   uLmais = u[i] + 0.5*Limiter(r_i)*(u[i+1]-u[i])
                   uRmais = u[i+1] - 0.5*Limiter(r_i_mais)*(u[i+2]-u[i+1])
                   uLmenos = u[i-1] + 0.5*Limiter(r_i_menos)*(u[i]-u[i-1])
                   uRmenos = u[i] - 0.5*Limiter(r_i)*(u[i+1]-u[i])
                 end

    # Fluxo na face i+1/2 (direita)
    F_mais = vel[i+1] * (uLmais * (sign(vel[i+1]) + 1) / 2 + uRmais * (-sign(vel[i+1]) + 1) / 2)

    # Fluxo na face i-1/2 (esquerda)
    F_menos = vel[i] * (uLmenos * (sign(vel[i]) + 1) / 2 + uRmenos * (-sign(vel[i]) + 1) / 2)    

        U[i] = fi0[i] - (Δt/Δx)*(F_mais - F_menos) + Δt*(c4[i] + c5[i] + c6[i] + c7[i])
             end
            end  =#
  end
  U[n] = U[n-1]
  flow.uδD[3:(n+1), :, 1] .= U[2:n]
  flow.uδD[end, :, 1] = flow.uδD[end-1, :, 1]
  return (flow.uδD[:, :, 1])
end
# #hL = MUSCL(hLflow, 0, 0, 0, 0, 0, 0, meshP, CFL, Δt, Δx, 0, θ, 0, 0, 0, 0, fi0, P, αL0, αG0, ρL0, ρG0, uL0, uG0, PØ, hL0) 