# ==============================================================================
# closures_stratified.jl
#
# Closures geométricos e de atrito para escoamento horizontal estratificado.
#
# Objetivo:
#   Substituir o bloco de closures vertical/anelar do solver quando:
#       inp.energyflag == 0
#       inp.flowpattern == 2   # sugestão: 2 = horizontal estratificado
#
# Este arquivo usa a relação alpha_L -> h_L já resolvida por tabela rápida.
# Pode ser usado junto com hLCompute.jl ou de forma independente, pois as funções
# de inversão alpha_frac -> hL também estão incluídas aqui.
#
# Convenção:
#   alpha_frac = fração volumétrica líquida do solver, alpha_L ∈ [0,1]
#   hL         = altura da lâmina líquida
#   D          = diâmetro do duto
#
# Observação importante:
#   O flow.jl atual usa um único vetor τw nas equações de momento. Neste arquivo,
#   calculamos τwL e τwG separadamente, mas mapeamos τw .= τwL para manter
#   compatibilidade imediata com o solver atual. Para incluir corretamente o
#   atrito gás-parede na equação do gás, o rhs! em flow.jl precisa ser ajustado.
# ==============================================================================

# ------------------------------------------------------------------------------
# Inversão rápida alpha_frac -> hL
# ------------------------------------------------------------------------------

@inline function alpha_geom_eta(eta::Float64)
    eta = clamp(eta, 0.0, 1.0)
    z = clamp(2.0 * eta - 1.0, -1.0, 1.0)
    return 0.5 * (π - acos(z) + z * sqrt(max(0.0, 1.0 - z*z)))
end

@inline function frac_to_geom(alpha_frac::Float64)
    return (π / 2) * clamp(alpha_frac, 0.0, 1.0)
end

struct FastAlphaFracTable
    eta::Vector{Float64}
    n::Int
    alpha_geom_min::Float64
    alpha_geom_max::Float64
    dalpha::Float64
    inv_dalpha::Float64
end

function build_fast_alpha_frac_table(n_alpha::Int=4001; n_aux::Int=20001)
    eta_aux = collect(range(0.0, 1.0, length=n_aux))
    alpha_aux = similar(eta_aux)

    @inbounds for i in eachindex(eta_aux)
        alpha_aux[i] = alpha_geom_eta(eta_aux[i])
    end

    alpha_geom_min = 0.0
    alpha_geom_max = π / 2
    dalpha = (alpha_geom_max - alpha_geom_min) / (n_alpha - 1)
    inv_dalpha = 1.0 / dalpha

    eta_table = Vector{Float64}(undef, n_alpha)
    j_aux = 2

    @inbounds for j in 1:n_alpha
        alpha_target = alpha_geom_min + (j - 1) * dalpha

        if alpha_target <= alpha_aux[1]
            eta_table[j] = eta_aux[1]
            continue
        elseif alpha_target >= alpha_aux[end]
            eta_table[j] = eta_aux[end]
            continue
        end

        while j_aux <= n_aux && alpha_aux[j_aux] < alpha_target
            j_aux += 1
        end

        j0 = j_aux - 1
        j1 = j_aux

        a0 = alpha_aux[j0]
        a1 = alpha_aux[j1]
        e0 = eta_aux[j0]
        e1 = eta_aux[j1]

        if abs(a1 - a0) < 1e-14
            eta_table[j] = 0.5 * (e0 + e1)
        else
            w = (alpha_target - a0) / (a1 - a0)
            eta_table[j] = e0 + w * (e1 - e0)
        end
    end

    return FastAlphaFracTable(
        eta_table,
        n_alpha,
        alpha_geom_min,
        alpha_geom_max,
        dalpha,
        inv_dalpha,
    )
end

@inline function invert_alpha_geom_fast(alpha_geom::Float64, table::FastAlphaFracTable)
    if alpha_geom <= table.alpha_geom_min
        return 0.0
    elseif alpha_geom >= table.alpha_geom_max
        return 1.0
    end

    x = (alpha_geom - table.alpha_geom_min) * table.inv_dalpha + 1.0
    i0 = Int(floor(x))
    i0 = clamp(i0, 1, table.n - 1)
    i1 = i0 + 1
    w = x - i0

    e0 = table.eta[i0]
    e1 = table.eta[i1]

    return e0 + w * (e1 - e0)
end

@inline function invert_alpha_frac_fast(alpha_frac::Float64, table::FastAlphaFracTable)
    return invert_alpha_geom_fast(frac_to_geom(alpha_frac), table)
end

@inline function compute_hL_from_alpha_frac_fast(alpha_frac::Float64, D::Float64, table::FastAlphaFracTable)
    return D * invert_alpha_frac_fast(alpha_frac, table)
end

# ------------------------------------------------------------------------------
# Estado auxiliar das closures estratificadas
# ------------------------------------------------------------------------------

mutable struct StratifiedClosureState
    hL::Vector{Float64}
    SL::Vector{Float64}
    SG::Vector{Float64}
    Si::Vector{Float64}
    DhL::Vector{Float64}
    DhG::Vector{Float64}
    dαdhL::Vector{Float64}
    ReL::Vector{Float64}
    ReG::Vector{Float64}
    fL::Vector{Float64}
    fG::Vector{Float64}
    fi::Vector{Float64}
    F::Vector{Float64}
    τwL::Vector{Float64}
    τwG::Vector{Float64}
    τi::Vector{Float64}
end

function allocate_stratified_closure_state(N::Int)
    return StratifiedClosureState(
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
        zeros(N),
    )
end

# ------------------------------------------------------------------------------
# Fórmulas geométricas do estratificado
# ------------------------------------------------------------------------------

@inline function stratified_geometry_scalar(αL::Float64, D::Float64, A::Float64, table::FastAlphaFracTable)
    αL_c = clamp(αL, 1e-10, 1.0 - 1e-10)
    αG_c = 1.0 - αL_c

    hL = compute_hL_from_alpha_frac_fast(αL_c, D, table)

    z = clamp(2.0 * hL / D - 1.0, -1.0, 1.0)
    root = sqrt(max(0.0, 1.0 - z*z))

    SL = D * (π - acos(z))
    SG = π * D - SL
    Si = D * root

    DhL = 4.0 * αL_c * A / max(SL, eps())
    DhG = 4.0 * αG_c * A / max(SG + Si, eps())

    # Da lâmina na PPTX:
    # dα_L/dh_L = (4/D^2) * sqrt(hL*(D - hL))
    dαdhL = (4.0 / D^2) * sqrt(max(0.0, hL * (D - hL)))

    return hL, SL, SG, Si, DhL, DhG, dαdhL
end

function compute_stratified_geometry!(
    C::StratifiedClosureState,
    αL::AbstractVector{Float64},
    D::Float64,
    A::Float64,
    table::FastAlphaFracTable,
)
    @assert length(C.hL) == length(αL)

    @inbounds for i in eachindex(αL)
        hL, SL, SG, Si, DhL, DhG, dαdhL = stratified_geometry_scalar(αL[i], D, A, table)
        C.hL[i] = hL
        C.SL[i] = SL
        C.SG[i] = SG
        C.Si[i] = Si
        C.DhL[i] = DhL
        C.DhG[i] = DhG
        C.dαdhL[i] = dαdhL
    end

    return nothing
end

# ------------------------------------------------------------------------------
# Atrito e cisalhamento
# ------------------------------------------------------------------------------

@inline function friction_factor_stratified(Re::Float64, Dh::Float64, roughness::Float64)
    Re_safe = max(abs(Re), 1e-12)
    Dh_safe = max(Dh, 1e-12)
    eps_rel = max(roughness / Dh_safe, 0.0)

    # f_k = max[16/Re, 0.001375{1 + [2e4 eps/Dh + 1e6/Re]^(1/3)}]
    lam = 16.0 / Re_safe
    turb = 0.001375 * (1.0 + (2.0e4 * eps_rel + 1.0e6 / Re_safe)^(1.0 / 3.0))

    return max(lam, turb)
end

#=
@inline function friction_factor_gas_for_interfacial(ReG::Float64)
    Re_safe = max(abs(ReG), 1e-12)
    return max(16.0 / Re_safe, 0.046 / (Re_safe^0.2))
end
=#

@inline function stratified_F(
    uG::Float64,
    ρL::Float64,
    ρG::Float64,
    dαdhL::Float64,
    αL::Float64,
    g::Float64,
    beta::Float64,
)
    denom = max((ρL - ρG) * max(1.0 - αL, 1e-10) * abs(g * cos(beta)), 1e-12)
    arg = (ρG / denom) * max(dαdhL, 0.0)
    return uG * sqrt(max(arg, 0.0))
end


@inline function interfacial_friction_factor(fG::Float64, F::Float64, hL::Float64, D::Float64)
    if F > 0.36
        return fG * (1.0 + 29.7 * (F - 0.36)^0.67 * (max(hL / D, 0.0))^0.2)
    else
        return fG
    end
end

function compute_stratified_friction!(
    C::StratifiedClosureState,
    ρL::AbstractVector{Float64},
    ρG::AbstractVector{Float64},
    uL::AbstractVector{Float64},
    uG::AbstractVector{Float64},
    αL::AbstractVector{Float64},
    μL::Float64,
    μG::Float64,
    D::Float64,
    g::Float64;
    roughness::Float64=0.0,
    beta::Float64=0.0,
)
    @assert length(C.hL) == length(ρL) == length(ρG) == length(uL) == length(uG) == length(αL)

    @inbounds for i in eachindex(ρL)
        ReL = ρL[i] * C.DhL[i] * abs(uL[i]) / max(μL, eps())
        ReG = ρG[i] * C.DhG[i] * abs(uG[i]) / max(μG, eps())

        fL = friction_factor_stratified(ReL, C.DhL[i], roughness)

        #=
        fG_wall = friction_factor_stratified(ReG, C.DhG[i], roughness)
        fG_int_base = friction_factor_gas_for_interfacial(ReG)

        F = stratified_F(uG[i], ρL[i], ρG[i], C.dαdhL[i], clamp(αL[i], 1e-10, 1.0 - 1e-10), g, beta)
        fi = interfacial_friction_factor(fG_int_base, F, C.hL[i], D)
        =#

        fG_wall = friction_factor_stratified(ReG, C.DhG[i], roughness)
        F = stratified_F(
            uG[i],
            ρL[i],
            ρG[i],
            C.dαdhL[i],
            clamp(αL[i], 1e-10, 1.0 - 1e-10),
            g,
            beta,
        )
        fi = interfacial_friction_factor(fG_wall, F, C.hL[i], D)

        C.ReL[i] = ReL
        C.ReG[i] = ReG
        C.fL[i] = fL
        C.fG[i] = fG_wall
        C.F[i] = F
        C.fi[i] = fi

        C.τwL[i] = 0.5 * fL * ρL[i] * uL[i] * abs(uL[i])
        C.τwG[i] = 0.5 * fG_wall * ρG[i] * uG[i] * abs(uG[i])
        C.τi[i]  = 0.5 * fi * ρG[i] * (uG[i] - uL[i]) * abs(uG[i] - uL[i])
    end

    return nothing
end

# ------------------------------------------------------------------------------
# Função principal para plugar no solver organizado
# ------------------------------------------------------------------------------

function compute_stratified_closures!(
    C::StratifiedClosureState,
    S,
    inp,
    grid,
    table::FastAlphaFracTable;
    roughness::Float64=0.0,
    beta::Float64=0.0,
)
    # Geometria a partir de alpha_L no centro das células
    compute_stratified_geometry!(C, S.αL, inp.D, grid.A, table)

    # Usa velocidades centradas. Chamar antes:
    #   S.uLC[i] = (S.uL[i] + S.uL[i+1])/2
    #   S.uGC[i] = (S.uG[i] + S.uG[i+1])/2
    compute_stratified_friction!(
        C,
        S.ρL,
        S.ρG,
        S.uLC,
        S.uGC,
        S.αL,
        inp.μL,
        inp.μG,
        inp.D,
        inp.g;
        roughness=roughness,
        beta=beta,
    )

    # Compatibilidade com arrays já existentes no solver
    S.ReL[:] .= C.ReL[:]
    S.ReG[:] .= C.ReG[:]
    S.fw[:]  .= C.fL[:]
    S.fi[:]  .= C.fi[:]

    # O solver atual usa τw único. Usamos τwL como aproximação compatível.
    # Para gás-parede, consultar C.τwG e adaptar rhs! em flow.jl se necessário.
    S.τwL[1:length(C.τwL)] .= C.τwL[:]
    S.τwG[1:length(C.τwG)] .= C.τwG[:]
    S.τi[1:length(C.τi)]  .= C.τi[:]

    # Condições de contorno simples para vetores em faces/células estendidas
    if length(S.τwL) == length(C.τwL) + 1
        S.τwL[end] = S.τwL[end-1]
    end

    if length(S.τwG) == length(C.τwG) + 1
        S.τwG[end] = S.τwG[end-1]
    end

    if length(S.τi) == length(C.τi) + 1
        S.τi[end] = S.τi[end-1]
    end

    if length(S.τi) >= 2
        S.τi[1] = S.τi[2]
    end

    return nothing
end

# ------------------------------------------------------------------------------
# Diagnóstico simples
# ------------------------------------------------------------------------------

function print_stratified_closure_summary(C::StratifiedClosureState; idx::Int=1)
    println("Stratified closure summary at cell ", idx)
    println("hL   = ", C.hL[idx])
    println("SL   = ", C.SL[idx])
    println("SG   = ", C.SG[idx])
    println("Si   = ", C.Si[idx])
    println("DhL  = ", C.DhL[idx])
    println("DhG  = ", C.DhG[idx])
    println("ReL  = ", C.ReL[idx])
    println("ReG  = ", C.ReG[idx])
    println("fL   = ", C.fL[idx])
    println("fG   = ", C.fG[idx])
    println("fi   = ", C.fi[idx])
    println("τwL  = ", C.τwL[idx])
    println("τwG  = ", C.τwG[idx])
    println("τi   = ", C.τi[idx])
    return nothing
end
