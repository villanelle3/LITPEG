function read_input(path::AbstractString="Entrada.txt")
    n_dados = 41
    dados = zeros(n_dados)
    texto = readlines(path)
    for i = 2:2:2*n_dados
        dados[Int(i ÷ 2)] = parse(Float64, texto[i])
    end
    return InputData(
        Comp = dados[1],
        D = dados[2],
        θ = dados[3],
        n_segmentos = Int(dados[4]),
        tempo = dados[5],
        uL_inlet = dados[6],
        uG_inlet = dados[7],
        ρLi = dados[8],
        ρGi = dados[9],
        P_out = dados[10],
        αG_inlet = dados[11],
        αG_0 = dados[12],
        T_intet = dados[13],
        T_0 = dados[14],
        CFL = dados[15],
        Tol_P = dados[16],
        Tol_u = dados[17],
        Tol_α = dados[18],
        Cp_L = dados[19],
        Cp_G = dados[20],
        kf_L = dados[21],
        kf_G = dados[22],
        T_ext = dados[23],
        μL = dados[24],
        μG = dados[25],
        kf_pipe = dados[26],
        D_ext = dados[27],
        g = dados[28],
        cG = dados[29],
        cL = dados[30],
        βG = dados[31],
        βL = dados[32],
        γ = dados[33],
        h_ext = dados[34],
        metodo = Int(dados[35]),
        Grau = Int(dados[36]),
        C = Int(dados[37]),
        energyflag = Int(dados[38]),
        flowpattern = Int(dados[39]),
        roughness = dados[40],
        estabmet = Int(dados[41])
    )
end

build_grid(inp::InputData) = GridData(
    N = inp.n_segmentos,
    Δx = inp.Comp / inp.n_segmentos,
    LMIN = 0.0,
    LMAX = inp.Comp,
    A = (π * inp.D^2) / 4,
    S = π * inp.D,
)

# ============================================================
# Alocação dos arrays
# ============================================================
