function allocate_state(inp::InputData, grid::GridData)
    N = grid.N
    # Tempo n+1
    P = zeros(N)
    αL = zeros(N)
    αG = zeros(N)
    αLf = zeros(N + 1)
    αGf = zeros(N + 1)
    uL = zeros(N + 1)
    uG = zeros(N + 1)
    ρL = zeros(N)
    ρG = zeros(N)
    ρLf = zeros(N + 1)
    ρGf = zeros(N + 1)
    hL = zeros(N)
    TL = zeros(N)
    hLf = zeros(N + 1)

    # Tempo n
    P0 = zeros(N)
    αL0 = zeros(N)
    αG0 = zeros(N)
    αLf0 = zeros(N + 1)
    αGf0 = zeros(N + 1)
    uL0 = zeros(N + 1)
    uG0 = zeros(N + 1)
    ρL0 = zeros(N)
    ρG0 = zeros(N)
    ρLf0 = zeros(N + 1)
    ρGf0 = zeros(N + 1)
    hL0 = zeros(N)
    hLf0 = zeros(N + 1)

    # Supostos / corretivos
    PØ = zeros(N)
    αLØ = zeros(N)
    αGØ = zeros(N)
    uLØ = zeros(N + 1)
    uGØ = zeros(N + 1)
    uLf = zeros(N)
    uGf = zeros(N)
    uLδ = zeros(N + 1)
    uGδ = zeros(N + 1)
    Pδ = zeros(N)
    Tδ = zeros(N)

    # Centros
    uLC = zeros(N)
    uGC = zeros(N)
    uLC0 = zeros(N)
    uGC0 = zeros(N)

    # Adimensionais
    ReL = zeros(N)
    ReG = zeros(N)
    R_núcleo = zeros(N)
    fw = zeros(N)
    fi = zeros(N)
    λ = zeros(N)
    τwL = zeros(N + 1)
    τwG = zeros(N + 1)
    τi = zeros(N + 1)
    Nu = zeros(N)
    ULM = zeros(N)
    U_conv = zeros(N)
    R_tot = zeros(N)

    fi0 = zeros(N)
    dif_art_L = zeros(N + 1)
    dif_art_G = zeros(N + 1)

    # Perímetros
    Si = zeros(N + 1)
    SwL = zeros(N + 1)
    SwG = zeros(N + 1)

    # Auxiliares
    Paux = zeros(N)
    αLaux = zeros(N)
    αGaux = zeros(N)
    uLaux = zeros(N + 1)
    uGaux = zeros(N + 1)

    a1 = zeros(N)
    a2 = zeros(N)
    a3 = zeros(N)
    a4 = zeros(N)
    a5 = zeros(N)
    a6 = zeros(N)
    a7 = zeros(N)
    a8 = zeros(N)

    # Pressão
    d_P = zeros(N)
    du_P = zeros(N)
    dl_P = zeros(N)
    b_P = zeros(N)

    return (
        P=P, αL=αL, αG=αG, αLf=αLf, αGf=αGf, uL=uL, uG=uG, ρL=ρL, ρG=ρG, ρLf=ρLf, ρGf=ρGf,
        hL=hL, TL=TL, hLf=hLf,
        P0=P0, αL0=αL0, αG0=αG0, αLf0=αLf0, αGf0=αGf0, uL0=uL0, uG0=uG0, ρL0=ρL0, ρG0=ρG0, ρLf0=ρLf0, ρGf0=ρGf0,
        hL0=hL0, hLf0=hLf0,
        PØ=PØ, αLØ=αLØ, αGØ=αGØ, uLØ=uLØ, uGØ=uGØ, uLf=uLf, uGf=uGf, uLδ=uLδ, uGδ=uGδ, Pδ=Pδ, Tδ=Tδ,
        uLC=uLC, uGC=uGC, uLC0=uLC0, uGC0=uGC0,
        ReL=ReL, ReG=ReG, R_núcleo=R_núcleo, fw=fw, fi=fi, λ=λ, τwL=τwL, τwG=τwG, τi=τi, SwL=SwL, SwG=SwG, Si=Si, Nu=Nu, 
        ULM=ULM, U_conv=U_conv, R_tot=R_tot,
        fi0=fi0, dif_art_L=dif_art_L, dif_art_G=dif_art_G,
        Paux=Paux, αLaux=αLaux, αGaux=αGaux, uLaux=uLaux, uGaux=uGaux,
        a1=a1, a2=a2, a3=a3, a4=a4, a5=a5, a6=a6, a7=a7, a8=a8,
        d_P=d_P, du_P=du_P, dl_P=dl_P, b_P=b_P,
    )
end

# ============================================================
# Inicialização global e do FR
# ============================================================
