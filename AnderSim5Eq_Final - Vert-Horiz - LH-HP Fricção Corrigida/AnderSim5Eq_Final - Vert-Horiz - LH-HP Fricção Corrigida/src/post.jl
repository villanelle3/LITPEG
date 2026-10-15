function finalize_results(S, inp::InputData, grid::GridData)
    df = DataFrame(
        Pr=S.P,
        UL=S.uL[1:grid.N],
        UG=S.uG[1:grid.N],
        AlfaL=S.αL,
        AlfaG=S.αG,
        Temperatura=S.TL,
        ρL=S.ρL,
        ρG=S.ρG,
    )
    CSV.write("export.csv", df)

    x = range(0, inp.Comp, length=grid.N)

    # =========================
    # Pressão
    # =========================
    pP = plot(
        x,
        S.P,
        label="P",
        xlabel="x (m)",
        ylabel="Pressão (Pa)",
        title="Pressão",
    )

    # =========================
    # Fração de gás
    # =========================
    pAlpha = plot(
        x,
        S.αG,
        label="αG",
        xlabel="x (m)",
        ylabel="Fração de gás",
        title="Fração volumétrica de gás",
    )

    # =========================
    # Velocidades
    # =========================
    pVel = plot(
        x,
        S.uL[1:grid.N],
        label="uL",
        xlabel="x (m)",
        ylabel="Velocidade (m/s)",
        title="Velocidades",
    )

    plot!(
        pVel,
        x,
        S.uG[1:grid.N],
        label="uG",
    )

    # =========================
    # Temperatura
    # =========================
    pT = plot(
        x,
        S.TL,
        label="T",
        xlabel="x (m)",
        ylabel="Temperatura (°C)",
        title="Temperatura",
    )

    # =========================
    # Layout final
    # =========================
    pFinal = plot(
        pP,
        pAlpha,
        pVel,
        pT,
        layout=(2, 2),
        size=(1200, 800),
    )

    display(pFinal)

    return df
end
