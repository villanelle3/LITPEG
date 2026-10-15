function initialize_globals!(inp::InputData, grid::GridData)
    global αG_inlet = inp.αG_inlet
    global αL_inlet = 1.0 - inp.αG_inlet
    global uL_inlet = inp.uL_inlet
    global uG_inlet = inp.uG_inlet
    global hL_intet = inp.Cp_L * inp.T_intet
    global hLi = inp.Cp_L * inp.T_0
    global Grau = inp.Grau
    global Λ = inp.Cp_L / inp.Cp_G
    global N = grid.N
    return nothing
end

function initialize_state!(S, inp::InputData, grid::GridData)
    αL_0 = 1.0 - inp.αG_0
    αL_in = 1.0 - inp.αG_inlet

    S.uLC[:] .= inp.uL_inlet
    S.uGC[:] .= inp.uG_inlet
    S.uLC0[:] .= inp.uL_inlet
    S.uGC0[:] .= inp.uG_inlet

    S.uL[:] .= inp.uL_inlet
    S.uG[:] .= inp.uG_inlet
    S.αL[:] .= αL_0
    S.αG[:] .= inp.αG_0
    S.ρL[:] .= inp.ρLi
    S.ρG[:] .= inp.ρGi
    S.P[:]  .= inp.P_out
    S.hL[:] .= hLi

    S.uL0[:] .= inp.uL_inlet
    S.uG0[:] .= inp.uG_inlet
    S.αL0[:] .= αL_in
    S.αG0[:] .= inp.αG_inlet
    S.ρL0[:] .= inp.ρLi
    S.ρG0[:] .= inp.ρGi
    S.P0[:]  .= inp.P_out
    S.hL0[:] .= hLi

    S.uLØ[:] .= inp.uL_inlet
    S.uGØ[:] .= inp.uG_inlet
    S.αLØ[:] .= S.αL0[:]
    S.αGØ[:] .= S.αG0[:]
    S.PØ[:]  .= inp.P_out

    S.d_P[1] = 1.0
    S.du_P[1] = 0.0
    S.dl_P[1] = 0.0
    S.b_P[1] = 0.0
    S.d_P[grid.N] = 1.0
    S.du_P[grid.N] = 0.0
    S.dl_P[grid.N] = 0.0
    S.b_P[grid.N] = 0.0
    return nothing
end

function initialize_fr!(inp::InputData, grid::GridData)
    meshP, innerPointsP = frgrid(grid.LMIN, grid.LMAX, grid.N, inp.Grau, inp.C)
    meshS, innerPointsS = frgrid(grid.LMIN - grid.Δx / 2, grid.LMAX + grid.Δx / 2, grid.N + 1, inp.Grau, inp.C)

    T = Float64

    f1(x) = inp.uL_inlet
    uL_IC = Function[f1]
    uLflow = Flow(grid.N + 1, inp.Grau, 1, 1; T=T)
    init!(uLflow, uL_IC, innerPointsS)
    uLflow.uδD[1, :, 1] .= inp.uL_inlet
    uLflow.uδD[end, :, 1] .= inp.uL_inlet
    uLδD = zeros(grid.N + 3, inp.Grau + 1)
    uLδD0 = zeros(grid.N + 3, inp.Grau + 1)
    uLδD[:, :] .= uLflow.uδD[:, :, 1]
    uLδD0[:, :] .= uLδD[:, :]

    f2(x) = inp.uG_inlet
    uG_IC = Function[f2]
    uGflow = Flow(grid.N + 1, inp.Grau, 1, 2; T=T)
    init!(uGflow, uG_IC, innerPointsS)
    uGflow.uδD[1, :, 1] .= inp.uG_inlet
    uGflow.uδD[end, :, 1] .= inp.uG_inlet
    uGδD = zeros(grid.N + 3, inp.Grau + 1)
    uGδD0 = zeros(grid.N + 3, inp.Grau + 1)
    uGδD[:, :] .= uGflow.uδD[:, :, 1]
    uGδD0[:, :] .= uGflow.uδD[:, :, 1]

    f3(x) = 1.0 - inp.αG_0
    αL_IC = Function[f3]
    αLflow = Flow(grid.N, inp.Grau, 2, 3; T=T)
    init!(αLflow, αL_IC, innerPointsP)
    αLflow.uδD[1, :, 1] .= 1.0 - inp.αG_inlet
    αLflow.uδD[end, :, 1] .= inp.αG_inlet
    αLδD = zeros(grid.N + 2, inp.Grau + 1)
    αLδD0 = zeros(grid.N + 2, inp.Grau + 1)
    αLδD[:, :] .= αLflow.uδD[:, :, 1]
    αLδD0[:, :] .= αLflow.uδD[:, :, 1]

    f4(x) = inp.αG_0
    αG_IC = Function[f4]
    αGflow = Flow(grid.N, inp.Grau, 2, 4; T=T)
    init!(αGflow, αG_IC, innerPointsP)
    αGflow.uδD[1, :, 1] .= inp.αG_inlet
    αGflow.uδD[end, :, 1] .= 1.0 - inp.αG_inlet
    αGδD = zeros(grid.N + 2, inp.Grau + 1)
    αGδD0 = zeros(grid.N + 2, inp.Grau + 1)
    αGδD[:, :] .= αGflow.uδD[:, :, 1]
    αGδD0[:, :] .= αGflow.uδD[:, :, 1]

    f5(x) = hLi
    hL_IC = Function[f5]
    hLflow = Flow(grid.N, inp.Grau, 3, 5; T=T)
    init!(hLflow, hL_IC, innerPointsP)
    hLflow.uδD[1, :, 1] .= hL_intet
    hLflow.uδD[end, :, 1] .= hL_intet
    hLδD = zeros(grid.N + 2, inp.Grau + 1)
    hLδD0 = zeros(grid.N + 2, inp.Grau + 1)
    hLδD[:, :] .= hLflow.uδD[:, :, 1]
    hLδD0[:, :] .= hLflow.uδD[:, :, 1]

    αLu = zeros(grid.N + 3, inp.Grau + 1)
    αGu = zeros(grid.N + 3, inp.Grau + 1)
    αLu0 = zeros(grid.N + 3, inp.Grau + 1)
    αGu0 = zeros(grid.N + 3, inp.Grau + 1)

    uLα = zeros(grid.N + 2, inp.Grau + 1)
    uGα = zeros(grid.N + 2, inp.Grau + 1)
    uLα0 = zeros(grid.N + 2, inp.Grau + 1)
    uGα0 = zeros(grid.N + 2, inp.Grau + 1)

    global V = vandermonde_legendre(inp.Grau, pontos)

    return (
        meshP=meshP, innerPointsP=innerPointsP, meshS=meshS, innerPointsS=innerPointsS,
        uLflow=uLflow, uGflow=uGflow, αLflow=αLflow, αGflow=αGflow, hLflow=hLflow,
        uLδD=uLδD, uLδD0=uLδD0, uGδD=uGδD, uGδD0=uGδD0,
        αLδD=αLδD, αLδD0=αLδD0, αGδD=αGδD, αGδD0=αGδD0,
        hLδD=hLδD, hLδD0=hLδD0,
        αLu=αLu, αGu=αGu, αLu0=αLu0, αGu0=αGu0,
        uLα=uLα, uGα=uGα, uLα0=uLα0, uGα0=uGα0
    )
end

# ============================================================
# Um passo temporal
# ============================================================
