function run_case(input_path::AbstractString="Entrada.txt"; verbose::Bool=true)
    inp = read_input(input_path) 
    grid = build_grid(inp)
    initialize_globals!(inp, grid)
    S = allocate_state(inp, grid)
    initialize_state!(S, inp, grid)
    FRD = initialize_fr!(inp, grid)
    alpha_table = build_fast_alpha_frac_table()
    strat = allocate_stratified_closure_state(grid.N)
    t = 0.0

    while t < inp.tempo
        aux1 = maximum(abs.(S.uL[:])) + inp.cL
        aux2 = maximum(abs.(S.uG[:])) + inp.cG
        wavevel = maximum([aux1, aux2])

        Δt = inp.metodo == 2 ? (inp.CFL * grid.Δx) / ((2 * inp.Grau + 1) * wavevel) :
                                      (inp.CFL * grid.Δx) / wavevel
        Δt = minimum([Δt, (inp.tempo - t)])
        t += Δt
        verbose && begin
            print(t)
            print("\n")
        end

        e_uL = 1.0
        e_uG = 1.0
        e_αL = 1.0
        e_αG = 1.0
        e_P = 1.0
        cont1 = 0

        while (e_uL > inp.Tol_u || e_uG > inp.Tol_u || e_αL > inp.Tol_α || e_αG > inp.Tol_α || e_P > inp.Tol_P) && cont1 < 1000
            cont1 += 1

            S.uLaux[:] = S.uLØ[:]
            S.uGaux[:] = S.uGØ[:]
            S.αLaux[:] = S.αLØ[:]
            S.αGaux[:] = S.αGØ[:]
            S.Paux[:] = S.PØ[:]

            # Faces
            S.αLf[1] = S.αL[1]; S.αGf[1] = S.αG[1]; S.ρLf[1] = S.ρL[1]; S.ρGf[1] = S.ρG[1]
            S.αLf0[1] = S.αL0[1]; S.αGf0[1] = S.αG0[1]; S.ρLf0[1] = S.ρL0[1]; S.ρGf0[1] = S.ρG0[1]
            S.αLf[grid.N+1] = S.αL[grid.N]; S.αGf[grid.N+1] = S.αG[grid.N]
            S.ρLf[grid.N+1] = S.ρL[grid.N]; S.ρGf[grid.N+1] = S.ρG[grid.N]
            S.αLf0[grid.N+1] = S.αL0[grid.N]; S.αGf0[grid.N+1] = S.αG0[grid.N]
            S.ρLf0[grid.N+1] = S.ρL0[grid.N]; S.ρGf0[grid.N+1] = S.ρG0[grid.N]

            for i = 2:grid.N
                S.αLf[i] = S.αL[i-1] * (sign(S.uLØ[i]) + 1) / 2 + S.αL[i] * (-sign(S.uLØ[i]) + 1) / 2
                S.αGf[i] = S.αG[i-1] * (sign(S.uGØ[i]) + 1) / 2 + S.αG[i] * (-sign(S.uGØ[i]) + 1) / 2
                S.αLf0[i] = S.αL0[i-1] * (sign(S.uL[i]) + 1) / 2 + S.αL0[i] * (-sign(S.uL[i]) + 1) / 2
                S.αGf0[i] = S.αG0[i-1] * (sign(S.uG[i]) + 1) / 2 + S.αG0[i] * (-sign(S.uG[i]) + 1) / 2
                S.ρLf[i] = S.ρL[i-1] * (sign(S.uL[i]) + 1) / 2 + S.ρL[i] * (-sign(S.uL[i]) + 1) / 2
                S.ρGf[i] = S.ρG[i-1] * (sign(S.uG[i]) + 1) / 2 + S.ρG[i] * (-sign(S.uG[i]) + 1) / 2
                S.ρLf0[i] = S.ρL0[i-1] * (sign(S.uL[i]) + 1) / 2 + S.ρL0[i] * (-sign(S.uL[i]) + 1) / 2
                S.ρGf0[i] = S.ρG0[i-1] * (sign(S.uG[i]) + 1) / 2 + S.ρG0[i] * (-sign(S.uG[i]) + 1) / 2
            end

            for i = 1:grid.N
                S.uLC[i] = (S.uL[i] + S.uL[i+1]) / 2
                S.uGC[i] = (S.uG[i] + S.uG[i+1]) / 2
                S.uLC0[i] = (S.uL0[i] + S.uL0[i+1]) / 2
                S.uGC0[i] = (S.uG0[i] + S.uG0[i+1]) / 2
            end

            FRD.αLu[1, :] .= S.αL[1]
            FRD.αGu[1, :] .= S.αG[1]
            FRD.αLu0[1, :] .= S.αL0[1]
            FRD.αGu0[1, :] .= S.αG0[1]
            FRD.αLu[end-1:end, :] .= S.αL[end]
            FRD.αGu[end-1:end, :] .= S.αG[end]
            FRD.αLu0[end-1:end, :] .= S.αL0[end]
            FRD.αGu0[end-1:end, :] .= S.αG0[end]

            for i = 2:grid.N+1
                for j = 1:inp.Grau+1
                    β1 = pontos[j] + 1
                    β2 = pontos[j] - 1
                    lj1 = lagrangeBasisAt(β1, pontos)
                    lj2 = lagrangeBasisAt(β2, pontos)
                    FRD.αLu[i, j] = sum(FRD.αLδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.αLδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                    FRD.αGu[i, j] = sum(FRD.αGδD[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.αGδD[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                    FRD.αLu0[i, j] = sum(FRD.αLδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.αLδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                    FRD.αGu0[i, j] = sum(FRD.αGδD0[i-1, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.αGδD0[i, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                end
            end

            if inp.flowpattern == 1 # Anular
                S.ReL[:] .= abs.((S.ρL[:] .* S.uLC[:] .* S.αL[:] * inp.D) ./ inp.μL)
                S.R_núcleo[:] .= (inp.D / 2) * (abs.((S.αG[:] .* S.uLC[:]) ./ (S.uGC[:] + S.αG[:] .* (S.uLC[:] - S.uGC[:])))) .^ (1 / 2)
                S.ReG[:] .= abs.(2 * (S.ρG[:] .* S.uGC[:] .* S.R_núcleo[:]) ./ inp.μG)

                S.fi[:] .= 0.079 ./ ((S.ReG[:] .+ 1e-12) .^ 0.25)
                S.τi[1:end-1] .= (S.fi[:] ./ 2) .* S.ρG[:] .* abs.(S.uGC[:] - S.uLC[:]) .* (S.uGC[:] - S.uLC[:])
                S.τi[1] = S.τi[2]

                S.fw[:] .= 0.046 ./ ((S.ReL[:] .+ 1e-12) .^ 0.2)
                S.λ[:] .= (1 .+ S.αG[:] + (2 * S.αG[:] .* log.(S.αG[:])) ./ (1 .- S.αG[:])) ./ (1 .- S.αG[:])
                S.τwL[1:end-1] .= S.fw[:] .* S.ρL[:] .* S.uLC[:] .* abs.(S.uLC[:]) ./ 2 .- (S.ρG[:] .- S.ρL[:]) .* (-inp.g) .* inp.D .* S.αG[:] .* S.λ[:] ./ 4
                S.τi[end] = S.τi[end-1]
                S.τwL[end] = S.τwL[end-1]

                S.Nu[:] .= ((S.fw[:] ./ 8) .* (S.ReL[:] .- 1000) .* ((inp.Cp_L * inp.μL) / inp.kf_L)) ./
                           (1.07 .+ 12.7 .* ((S.fw[:] ./ 8) .^ 0.5) .* ((((inp.Cp_L * inp.μL) / inp.kf_L) .^ (2 / 3)) .- 1))
                S.ULM[:] .= (inp.kf_L .* S.Nu[:]) ./ inp.D
            elseif inp.flowpattern == 2 # Estratificado
                compute_stratified_closures!(
                    strat,
                    S,
                    inp,
                    grid,
                    alpha_table;
                    roughness=inp.roughness,
                    beta=0.0,
                )
                S.τwL[1:grid.N] .= strat.τwL
                S.τwG[1:grid.N] .= strat.τwG
                S.τi[1:grid.N] .= strat.τi
                S.SwL[1:grid.N] .= strat.SL
                S.SwG[1:grid.N] .= strat.SG
                S.Si[1:grid.N] .= strat.Si
            end

            if inp.estabmet == 1
                for i = 2:grid.N
                    Δp_art = (inp.γ * S.αLf0[i] * S.ρLf0[i] * S.αGf0[i] * S.ρGf0[i] * (S.uG0[i] - S.uL0[i])^2) /
                             (S.αGf0[i] * S.ρLf0[i] + S.αLf0[i] * S.ρGf0[i])
                    S.dif_art_L[i] = Δp_art * (S.αL0[i] - S.αL0[i-1]) / grid.Δx
                    S.dif_art_G[i] = Δp_art * (S.αG0[i] - S.αG0[i-1]) / grid.Δx

                end
            elseif inp.estabmet == 2
                for i = 2:grid.N
                    S.dif_art_L[i] = ((S.P[i] + S.P[i-1]) / 2 - S.ρLf[i] * S.αLf[i] * cos(inp.θ) * (strat.dαdhL[i] + strat.dαdhL[i-1]) / 2) * (S.αL[i] - S.αL[i-1]) / grid.Δx
                    S.dif_art_G[i] = ((S.P[i] + S.P[i-1]) / 2) * (S.αG[i] - S.αG[i-1]) / grid.Δx
                end
            end

            if inp.metodo == 1
                FRD.uLδD[:, :] .= MUSCL(FRD.uLflow, FRD.uLδD0, S.αLf, S.ρLf, S.αLf0, S.ρLf0, S.uL, FRD.meshS, inp.CFL, Δt, grid.Δx, S.P, inp.θ, S.dif_art_L, S.uLC0, S.αL, S.ρL, S.uL0, S.τi, S.τwL, S.αGf0, inp.D, inp.flowpattern, S.SwL, S.Si)
                FRD.uGδD[:, :] .= MUSCL(FRD.uGflow, FRD.uGδD0, S.αGf, S.ρGf, S.αGf0, S.ρGf0, S.uG, FRD.meshS, inp.CFL, Δt, grid.Δx, S.P, inp.θ, S.dif_art_G, S.uGC0, S.αG, S.ρG, S.uG0, S.τi, S.τwG, S.αGf0, inp.D, inp.flowpattern, S.SwG, S.Si)
            else
                FRD.uLδD[:, :] .= FR(FRD.uLflow, FRD.uLδD0, FRD.αLu, S.ρLf, FRD.αLu0, S.ρLf0, FRD.uLα, FRD.meshS, inp.CFL, Δt, grid.Δx, S.P, inp.θ, S.dif_art_L, S.uLC, S.αL, S.ρL, S.uL0, S.τi, S.τwL, S.αGf0, inp.D, S.a1, S.a2, S.a3, S.a4, S.a5, S.a6, S.a7, S.a8, S.αG, S.ρG, S.uGC, S.αL, S.ρL, S.uLC, wavevel, inp.flowpattern, S.SwL, S.Si)
                FRD.uGδD[:, :] .= FR(FRD.uGflow, FRD.uGδD0, FRD.αGu, S.ρGf, FRD.αGu0, S.ρGf0, FRD.uGα, FRD.meshS, inp.CFL, Δt, grid.Δx, S.P, inp.θ, S.dif_art_G, S.uGC, S.αG, S.ρG, S.uG0, S.τi, S.τwG, S.αGf0, inp.D, S.a1, S.a2, S.a3, S.a4, S.a5, S.a6, S.a7, S.a8, S.αG, S.ρG, S.uGC, S.αL, S.ρL, S.uLC, wavevel, inp.flowpattern, S.SwG, S.Si)
            end

            FRD.uLδD[3, :] .= (FRD.uLδD[2, :] + FRD.uLδD[4, :]) ./ 2
            FRD.uGδD[3, :] .= (FRD.uGδD[2, :] + FRD.uGδD[4, :]) ./ 2
            FRD.uLflow.uδD[:, :, 1] .= FRD.uLδD[:, :]
            FRD.uGflow.uδD[:, :, 1] .= FRD.uGδD[:, :]

            for i = 1:grid.N+1
                S.uLØ[i] = sum(FRD.uLδD[i+1, :] .* weights[:]) / 2
                S.uGØ[i] = sum(FRD.uGδD[i+1, :] .* weights[:]) / 2
            end

            # Pressão
            S.d_P[grid.N] = 1.0
            S.du_P[grid.N] = 0.0
            S.dl_P[grid.N] = 0.0
            S.b_P[grid.N] = 0.0

            for i = 2:(grid.N-1)
                if S.uLØ[i+1] >= 0
                    αLen = S.αL0[i]; ρLen = S.ρL0[i]
                    αLwn = S.αL0[i-1]; ρLwn = S.ρL0[i-1]
                    αLe = S.αLØ[i]; ρLe = S.ρL[i]
                    αLw = S.αLØ[i-1]; ρLw = S.ρL[i-1]
                    a_Le = αLen * ρLen / Δt + αLen * ρLen * S.uL0[i+1] / grid.Δx
                    a_Lw = αLwn * ρLwn / Δt + αLwn * ρLwn * S.uL0[i] / grid.Δx
                else
                    αLen = S.αL0[i+1]; ρLen = S.ρL0[i+1]
                    αLwn = S.αL0[i]; ρLwn = S.ρL0[i]
                    αLe = S.αLØ[i+1]; ρLe = S.ρL[i+1]
                    αLw = S.αLØ[i]; ρLw = S.ρL[i]
                    a_Le = αLen * ρLen / Δt - αLen * ρLen * S.uL0[i+1] / grid.Δx
                    a_Lw = αLwn * ρLwn / Δt - αLwn * ρLwn * S.uL0[i] / grid.Δx
                end

                if S.uGØ[i+1] >= 0
                    αGen = S.αG0[i]; ρGen = S.ρG0[i]
                    αGwn = S.αG0[i-1]; ρGwn = S.ρG0[i-1]
                    αGe = S.αGØ[i]; ρGe = S.ρG[i]
                    αGw = S.αGØ[i-1]; ρGw = S.ρG[i-1]
                    a_Ge = αGen * ρGen / Δt + αGen * ρGen * S.uG0[i+1] / grid.Δx
                    a_Gw = αGwn * ρGwn / Δt + αGwn * ρGwn * S.uG0[i] / grid.Δx
                else
                    αGen = S.αG0[i+1]; ρGen = S.ρG0[i+1]
                    αGwn = S.αG0[i]; ρGwn = S.ρG0[i]
                    αGe = S.αGØ[i+1]; ρGe = S.ρG[i+1]
                    αGw = S.αGØ[i]; ρGw = S.ρG[i]
                    a_Ge = αGen * ρGen / Δt - αGen * ρGen * S.uG0[i+1] / grid.Δx
                    a_Gw = αGwn * ρGwn / Δt - αGwn * ρGwn * S.uG0[i] / grid.Δx
                end

                mL = (αLe * ρLe * S.uLØ[i+1] - αLw * ρLw * S.uLØ[i]) / inp.ρLi
                mG = (αGe * ρGe * S.uGØ[i+1] - αGw * ρGw * S.uGØ[i]) / inp.ρGi
                coef = (grid.Δx / Δt) * (S.αLØ[i] / (inp.ρLi * inp.cL^2) + S.αGØ[i] / (inp.ρGi * inp.cG^2))
                S.d_P[i] = coef + αLe * ρLe * (αLe / (grid.Δx * a_Le)) / inp.ρLi + αGe * ρGe * (αGe / (grid.Δx * a_Ge)) / inp.ρGi + (αLw * ρLw * (αLw / (grid.Δx * a_Lw)) / inp.ρLi + αGw * ρGw * (αGw / (grid.Δx * a_Gw)) / inp.ρGi)
                S.du_P[i] = -(αLe * ρLe * (αLe / (grid.Δx * a_Le)) / inp.ρLi + αGe * ρGe * (αGe / (grid.Δx * a_Ge)) / inp.ρGi)
                S.dl_P[i] = -(αLw * ρLw * (αLw / (grid.Δx * a_Lw)) / inp.ρLi + αGw * ρGw * (αGw / (grid.Δx * a_Gw)) / inp.ρGi)
                S.b_P[i] = -(mL + mG + (grid.Δx / Δt) * ((S.αLØ[i] * S.ρL[i] - S.αL0[i] * S.ρL0[i]) / inp.ρLi + (S.αGØ[i] * S.ρG[i] - S.αG0[i] * S.ρG0[i]) / inp.ρGi))
            end

            S.Pδ[:] = TDMA(S.d_P, S.du_P, S.dl_P, S.b_P)
            #S.Pδ[1] = S.Pδ[2]
            S.Pδ[1] = 2*S.Pδ[2] - S.Pδ[3]

            for i = 1:(grid.N-1)
                if S.uLØ[i+1] >= 0
                    αLe = S.αL0[i]; ρLe = S.ρL0[i]
                    a_L = αLe * ρLe / Δt + αLe * ρLe * S.uL0[i+1] / grid.Δx
                else
                    αLe = S.αL0[i+1]; ρLe = S.ρL0[i+1]
                    a_L = αLe * ρLe / Δt - αLe * ρLe * S.uL0[i+1] / grid.Δx
                end

                if S.uGØ[i+1] >= 0
                    αGe = S.αG0[i]; ρGe = S.ρG0[i]
                    a_G = αGe * ρGe / Δt + αGe * ρGe * S.uG0[i+1] / grid.Δx
                else
                    αGe = S.αG0[i+1]; ρGe = S.ρG0[i+1]
                    a_G = αGe * ρGe / Δt - αGe * ρGe * S.uG0[i+1] / grid.Δx
                end

                de_L = αLe / (grid.Δx * a_L)
                de_G = αGe / (grid.Δx * a_G)
                S.uLδ[i+1] = de_L * (S.Pδ[i] - S.Pδ[i+1])
                S.uGδ[i+1] = de_G * (S.Pδ[i] - S.Pδ[i+1])
            end

            S.P[:] = S.PØ[:] + S.Pδ[:]
            S.uL .= S.uLØ + S.uLδ
            S.uG .= S.uGØ + S.uGδ

            for i = 2:(grid.N+2)
                FRD.uLflow.uδD[i, :, 1] .+= S.uLδ[i-1]
                FRD.uGflow.uδD[i, :, 1] .+= S.uGδ[i-1]
            end
            FRD.uLflow.uδD[grid.N+3, :, 1] .= FRD.uLflow.uδD[grid.N+2, :, 1]
            FRD.uGflow.uδD[grid.N+3, :, 1] .= FRD.uGflow.uδD[grid.N+2, :, 1]

            S.uL[grid.N+1] = S.uL[grid.N]
            S.uG[grid.N+1] = S.uG[grid.N]

            for i = 2:grid.N
                S.Tδ[i] = (1 / inp.Cp_L) * (S.hL[i] - S.hL0[i])
            end

            for i = 1:grid.N
                #Líquido incompressível:# S.ρL[i] = S.ρL[i] + (1 / inp.cL^2) * S.Pδ[i] - S.ρL[i] * inp.βL * S.Tδ[i]
                #S.ρG[i] = S.ρG[i] + (1 / inp.cG^2) * S.Pδ[i] - S.ρG[i] * inp.βG * S.Tδ[i]
                S.ρG[i] = S.P[i]/(420*293)
            end

            # Reconstrução uLα / uGα
            for i = 1:grid.N+2
                for j = 1:inp.Grau+1
                    β1 = pontos[j] + 1
                    β2 = pontos[j] - 1
                    lj1 = lagrangeBasisAt(β1, pontos)
                    lj2 = lagrangeBasisAt(β2, pontos)
                    FRD.uLα[i, j] = sum(FRD.uLflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.uLflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                    FRD.uGα[i, j] = sum(FRD.uGflow.uδD[i, :, 1] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.uGflow.uδD[i+1, :, 1] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                    FRD.uLα0[i, j] = sum(FRD.uLδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.uLδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                    FRD.uGα0[i, j] = sum(FRD.uGδD0[i, :] .* lj1[:]) * 0.5 * (1 - sign(pontos[j])) + sum(FRD.uGδD0[i+1, :] .* lj2[:]) * 0.5 * (1 + sign(pontos[j]))
                end
            end

            for i = 1:grid.N
                S.uLC[i] = (S.uL[i] + S.uL[i+1]) / 2
                S.uGC[i] = (S.uG[i] + S.uG[i+1]) / 2
                S.uLC0[i] = (S.uL0[i] + S.uL0[i+1]) / 2
                S.uGC0[i] = (S.uG0[i] + S.uG0[i+1]) / 2
            end

            if inp.metodo == 1
                FRD.αLδD[:, :] .= MUSCL(FRD.αLflow, FRD.αLδD0, S.uLC, S.ρL0, S.uLC0, S.ρL0, S.uL, FRD.meshP, inp.CFL, Δt, grid.Δx, S.P, inp.θ, 0, S.αLf, S.αL, S.ρL, S.αL0, S.τi, S.τwL, S.αGf0, inp.D, inp.flowpattern, S.SwL, S.Si)
                #FRD.αGδD[:, :] .= MUSCL(FRD.αGflow, FRD.αGδD0, S.uGC, S.ρG0, S.uGC0, S.ρG0, S.uG, FRD.meshP, inp.CFL, Δt, grid.Δx, S.P, inp.θ, 0, S.αGf, S.αG, S.ρG, S.αG0, S.τi, S.τwG, S.αGf0, inp.D, inp.flowpattern, S.SwG, S.Si)
            else
                FRD.αLδD[:, :] .= FR(FRD.αLflow, FRD.αLδD0, FRD.uLα, S.ρL0, FRD.uLα0, S.ρL0, FRD.uLα, FRD.meshP, inp.CFL, Δt, grid.Δx, S.P, inp.θ, 0, S.αLf, S.αL, S.ρL, S.αL0, S.τi, S.τwL, S.αGf0, inp.D, S.a1, S.a2, S.a3, S.a4, S.a5, S.a6, S.a7, S.a8, S.αG, S.ρG, S.uGC, S.αL, S.ρL, S.uLC, wavevel, inp.flowpattern, S.SwL, S.Si)
                #FRD.αGδD[:, :] .= FR(FRD.αGflow, FRD.αGδD0, FRD.uGα, S.ρG0, FRD.uGα0, S.ρG0, FRD.uGα, FRD.meshP, inp.CFL, Δt, grid.Δx, S.P, inp.θ, 0, S.αGf, S.αG, S.ρG, S.αG0, S.τi, S.τwG, S.αGf0, inp.D, S.a1, S.a2, S.a3, S.a4, S.a5, S.a6, S.a7, S.a8, S.αG, S.ρG, S.uGC, S.αL, S.ρL, S.uLC, wavevel, inp.flowpattern, S.SwG, S.Si)
            end

            for j = 2:grid.N
                for k = 1:(inp.Grau+1)
                    #FRD.αLδD[j, k] = FRD.αLδD[j, k] / (FRD.αLδD[j, k] + FRD.αGδD[j, k])
                    FRD.αGδD[j, k] = 1.0 - FRD.αLδD[j, k]
                end
            end

            for i = 1:grid.N
                S.αL[i] = sum(FRD.αLδD[i+1, :] .* weights[:]) / 2
                S.αG[i] = sum(FRD.αGδD[i+1, :] .* weights[:]) / 2
            end
                S.αL[end] = 2*S.αL[end-1] - S.αL[end-2]
                S.αG[end] = 2*S.αG[end-1] - S.αG[end-2]

            # Faces com h
            S.αLf[1] = S.αL[1]; S.αGf[1] = S.αG[1]; S.ρLf[1] = S.ρL[1]; S.ρGf[1] = S.ρG[1]
            S.αLf0[1] = S.αL0[1]; S.αGf0[1] = S.αG0[1]; S.ρLf0[1] = S.ρL0[1]; S.ρGf0[1] = S.ρG0[1]
            S.hLf[1] = S.hL[1]; S.hLf0[1] = S.hL0[1]
            S.αLf[grid.N+1] = S.αL[grid.N]; S.αGf[grid.N+1] = S.αG[grid.N]
            S.ρLf[grid.N+1] = S.ρL[grid.N]; S.ρGf[grid.N+1] = S.ρG[grid.N]
            S.αLf0[grid.N+1] = S.αL0[grid.N]; S.αGf0[grid.N+1] = S.αG0[grid.N]
            S.ρLf0[grid.N+1] = S.ρL0[grid.N]; S.ρGf0[grid.N+1] = S.ρG0[grid.N]
            S.hLf[grid.N+1] = S.hL[grid.N]; S.hLf0[grid.N+1] = S.hL0[grid.N]

            for i = 2:grid.N
                S.αLf[i] = S.αL[i-1] * (sign(S.uLØ[i]) + 1) / 2 + S.αL[i] * (-sign(S.uLØ[i]) + 1) / 2
                S.αGf[i] = S.αG[i-1] * (sign(S.uGØ[i]) + 1) / 2 + S.αG[i] * (-sign(S.uGØ[i]) + 1) / 2
                S.αLf0[i] = S.αL0[i-1] * (sign(S.uL[i]) + 1) / 2 + S.αL0[i] * (-sign(S.uL[i]) + 1) / 2
                S.αGf0[i] = S.αG0[i-1] * (sign(S.uG[i]) + 1) / 2 + S.αG0[i] * (-sign(S.uG[i]) + 1) / 2
                S.ρLf[i] = S.ρL[i-1] * (sign(S.uL[i]) + 1) / 2 + S.ρL[i] * (-sign(S.uL[i]) + 1) / 2
                S.ρGf[i] = S.ρG[i-1] * (sign(S.uG[i]) + 1) / 2 + S.ρG[i] * (-sign(S.uG[i]) + 1) / 2
                S.ρLf0[i] = S.ρL0[i-1] * (sign(S.uL[i]) + 1) / 2 + S.ρL0[i] * (-sign(S.uL[i]) + 1) / 2
                S.ρGf0[i] = S.ρG0[i-1] * (sign(S.uG[i]) + 1) / 2 + S.ρG0[i] * (-sign(S.uG[i]) + 1) / 2
                S.hLf[i] = S.hL[i-1] * (sign(S.uL[i]) + 1) / 2 + S.hL[i] * (-sign(S.uL[i]) + 1) / 2
                S.hLf0[i] = S.hL0[i-1] * (sign(S.uL[i]) + 1) / 2 + S.hL0[i] * (-sign(S.uL[i]) + 1) / 2
            end
            
            
            # Energia
        if inp.energyflag == 1
            S.U_conv[:] .= S.ULM[:] .* (1 ./ (1 .- S.αG[:])) .^ 0.9
            S.R_tot[:] .= (1 ./ (S.U_conv[:] .* inp.D * π .* grid.Δx)) .+ log(inp.D_ext / inp.D) ./ (2 * π * inp.kf_pipe * grid.Δx) .+ (1 ./ (inp.h_ext .* (grid.Δx .* π .* inp.D_ext)))

            u = zeros(grid.N)
            u[:] = S.hL[:]
            S.hL[1] = hL_intet

            for i = 2:grid.N
                if i == 2
                    r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
                    r_i_menos = 0
                    r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
                    uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                    uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
                    uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
                    uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                elseif i == (grid.N - 1)
                    r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
                    r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
                    r_i_mais = 0
                    uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                    uRmais = u[i+1]
                    uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
                    uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                elseif i == grid.N
                    r_i = (u[i] - u[i-1]) / (u[i] - u[i] + eps())
                    r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
                    r_i_mais = 0
                    uLmais = u[i]
                    uRmais = u[i]
                    uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
                    uRmenos = u[i]
                else
                    r_i = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps())
                    r_i_menos = (u[i-1] - u[i-2]) / (u[i] - u[i-1] + eps())
                    r_i_mais = (u[i+1] - u[i]) / (u[i+2] - u[i+1] + eps())
                    uLmais = u[i] + 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                    uRmais = u[i+1] - 0.5 * Limiter(r_i_mais) * (u[i+2] - u[i+1])
                    uLmenos = u[i-1] + 0.5 * Limiter(r_i_menos) * (u[i] - u[i-1])
                    uRmenos = u[i] - 0.5 * Limiter(r_i) * (u[i+1] - u[i])
                end

                F_mais = ((S.αLf[i+1] * S.ρLf[i+1] * S.uL[i+1] * (uLmais + 0.5 * S.uL[i+1]^2))+(S.αGf[i+1] * S.ρGf[i+1] * S.uG[i+1] * (Λ * uLmais + 0.5 * S.uG[i+1]^2)))*(sign((S.uL[i]+S.uL[i+1])/2)+1)/2 +
                         ((S.αLf[i] * S.ρLf[i] * S.uL[i] * (uRmais + 0.5 * S.uL[i]^2))+(S.αGf[i] * S.ρGf[i] * S.uG[i] * (Λ * uRmais + 0.5 * S.uG[i]^2)))*(-sign((S.uL[i]+S.uL[i+1])/2)+1)/2

                F_menos = ((S.αLf[i+1] * S.ρLf[i+1] * S.uL[i+1] * (uLmenos + 0.5 * S.uL[i+1]^2))+(S.αGf[i+1] * S.ρGf[i+1] * S.uG[i+1] * (Λ * uLmenos + 0.5 * S.uG[i+1]^2)))*(sign((S.uL[i]+S.uL[i+1])/2)+1)/2 +
                          ((S.αLf[i] * S.ρLf[i] * S.uL[i] * (uRmenos + 0.5 * S.uL[i]^2))+(S.αGf[i] * S.ρGf[i] * S.uG[i] * (Λ * uRmenos + 0.5 * S.uG[i]^2)))*(-sign((S.uL[i]+S.uL[i+1])/2)+1)/2

                S.a3[i] = (-Δt / grid.Δx) * (F_mais - F_menos)
                S.a1[i] = -0.5 * (S.αG[i] * S.ρG[i] * S.uG[i]^2 + S.αL[i] * S.ρL[i] * S.uL[i]^2)
                S.a2[i] = S.αG0[i] * S.ρG0[i] * (Λ * S.hL0[i] + 0.5 * S.uG0[i]^2) + S.αL0[i] * S.ρL0[i] * (S.hL0[i] + 0.5 * S.uL0[i]^2)
                S.a4[i] = S.P[i] - S.PØ[i]
                S.a5[i] = Δt * inp.g * sin(inp.θ) * (S.αL0[i] * S.ρL0[i] * S.uL0[i] + S.αG0[i] * S.ρG0[i] * S.uG0[i])
                an = S.hL0[i] / inp.Cp_L
                if an >= inp.T_ext
                    S.a6[i] = (Δt * S.hL0[i] * grid.S) / (inp.Cp_L * S.R_tot[i] * grid.A)
                    S.a7[i] = (-inp.T_ext * Δt * grid.S) / (S.R_tot[i] * grid.A)
                else
                    S.a6[i] = (-Δt * S.hL0[i] * grid.S) / (inp.Cp_L * S.R_tot[i] * grid.A)
                    S.a7[i] = (inp.T_ext * Δt * grid.S) / (S.R_tot[i] * grid.A)
                end
                S.a8[i] = S.αG[i] * S.ρG[i] * Λ + S.αL[i] * S.ρL[i]
            end

            if inp.metodo == 1
                S.hL[2:end] = (S.a1[2:end] + S.a2[2:end] + S.a3[2:end] + S.a4[2:end] + S.a5[2:end] + S.a6[2:end] + S.a7[2:end]) ./ S.a8[2:end]
            else
                FRD.hLδD[:, :] .= FR(FRD.hLflow, FRD.hLδD0, S.uLC, S.ρL0, S.uLC0, S.ρL0, S.uL, FRD.meshP, inp.CFL, Δt, grid.Δx, S.P, inp.θ, 0, S.αLf, S.αL, S.ρL, S.αL0, S.τi, S.τwL, S.αGf0, inp.D, S.a1, S.a2, S.a3, S.a4, S.a5, S.a6, S.a7, S.a8, S.αG, S.ρG, S.uGC, S.αL, S.ρL, S.uLC, wavevel, inp.flowpattern, S.SwL, S.Si)
                FRD.hLflow.uδD[:, :, 1] .= FRD.hLδD[:, :]
                for i = 2:grid.N
                    S.hL[i] = sum(FRD.hLδD[i+1, :] .* weights[:]) / 2
                end
            end
        else
            S.hL .= S.hL0
            S.TL .= S.TL
        end

            e_uL = norm(S.uL - S.uLaux) / norm(S.uLaux)
            e_uG = norm(S.uG - S.uGaux) / norm(S.uGaux)
            e_αL = norm(S.αL - S.αLaux) / norm(S.αLaux)
            e_αG = norm(S.αG - S.αGaux) / norm(S.αGaux)
            e_P = norm(S.P - S.Paux) / norm(S.Paux)

            S.uLØ[:] = S.uL[:]
            S.uGØ[:] = S.uG[:]
            S.αLØ[:] = S.αL[:]
            S.αGØ[:] = S.αG[:]
            S.PØ[:] = S.P[:]
        end

        S.uL0[:] = S.uL[:]
        S.uG0[:] = S.uG[:]
        S.αL0[:] = S.αL[:]
        S.αG0[:] = S.αG[:]
        S.P0[:] = S.P[:]
        S.ρG0[:] = S.ρG[:]
        S.ρL0[:] = S.ρL[:]
        S.hL0[:] = S.hL[:]
        S.TL[:] = S.hL[:] / inp.Cp_L

        FRD.uLδD0[:, :] .= FRD.uLδD[:, :]
        FRD.uGδD0[:, :] .= FRD.uGδD[:, :]
        FRD.αLδD0[:, :] .= FRD.αLδD[:, :]
        FRD.αGδD0[:, :] .= FRD.αGδD[:, :]
    end

    return finalize_results(S, inp, grid)
end

