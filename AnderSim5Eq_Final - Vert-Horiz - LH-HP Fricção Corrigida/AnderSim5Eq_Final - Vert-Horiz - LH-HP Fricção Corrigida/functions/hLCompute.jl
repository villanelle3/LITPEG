# ============================================================
# Inversão rápida alpha_frac -> h_L
# usando apenas lookup table + interpolação linear
#
# Entrada do usuário:
#   alpha_frac ∈ [0, 1]
#
# Conversão interna:
#   alpha_geom = (π/2) * alpha_frac
#
# Geometria usada:
#   alpha_geom = (1/2) * [ π - acos(2hL/D - 1)
#                        + (2hL/D - 1)*sqrt(1 - (2hL/D - 1)^2) ]
# ============================================================

# ------------------------------------------------------------
# Alfa geométrico em função de eta = hL/D
# Faixa:
#   eta = 0   -> alpha_geom = 0
#   eta = 1   -> alpha_geom = π/2
# ------------------------------------------------------------
@inline function alpha_geom_eta(eta::Float64)
    eta = clamp(eta, 0.0, 1.0)
    z = clamp(2.0*eta - 1.0, -1.0, 1.0)
    return 0.5 * (π - acos(z) + z * sqrt(max(0.0, 1.0 - z*z)))
end

# ------------------------------------------------------------
# Conversão da fração líquida do simulador para alfa geométrico
# alpha_frac ∈ [0,1]  -> alpha_geom ∈ [0, π/2]
# ------------------------------------------------------------
@inline function frac_to_geom(alpha_frac::Float64)
    return (π/2) * clamp(alpha_frac, 0.0, 1.0)
end

# ============================================================
# Estrutura da tabela rápida
# Armazena eta(alpha_geom) em uma malha uniforme de alpha_geom
# ============================================================
struct FastAlphaFracTable
    eta::Vector{Float64}
    n::Int
    alpha_geom_min::Float64
    alpha_geom_max::Float64
    dalpha::Float64
    inv_dalpha::Float64
end

# ============================================================
# Construção da tabela
#
# Estratégia:
# 1) cria uma tabela auxiliar fina em eta
# 2) calcula alpha_geom(eta)
# 3) monta uma tabela uniforme em alpha_geom
# 4) interpola eta para cada alpha_geom uniforme
# ============================================================
function build_fast_alpha_frac_table(n_alpha::Int=4001; n_aux::Int=20001)

    # Tabela auxiliar uniforme em eta
    eta_aux = collect(range(0.0, 1.0, length=n_aux))
    alpha_aux = similar(eta_aux)

    @inbounds for i in eachindex(eta_aux)
        alpha_aux[i] = alpha_geom_eta(eta_aux[i])
    end

    alpha_geom_min = 0.0
    alpha_geom_max = π/2
    dalpha = (alpha_geom_max - alpha_geom_min) / (n_alpha - 1)
    inv_dalpha = 1.0 / dalpha

    eta_table = Vector{Float64}(undef, n_alpha)

    # Monta tabela eta(alpha_geom)
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
        inv_dalpha
    )
end

# ============================================================
# Inversão rápida: alpha_geom -> eta
# Usa índice direto + interpolação linear
# ============================================================
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

# ============================================================
# Função principal:
# recebe alpha_frac ∈ [0,1] e retorna eta = hL/D
# ============================================================
@inline function invert_alpha_frac_fast(alpha_frac::Float64, table::FastAlphaFracTable)
    alpha_geom = frac_to_geom(alpha_frac)
    return invert_alpha_geom_fast(alpha_geom, table)
end

# ============================================================
# Função principal:
# recebe alpha_frac ∈ [0,1] e retorna hL
# ============================================================
@inline function compute_hL_from_alpha_frac_fast(
    alpha_frac::Float64,
    D::Float64,
    table::FastAlphaFracTable
)
    eta = invert_alpha_frac_fast(alpha_frac, table)
    return D * eta
end

# ============================================================
# Atualização vetorial - D constante
# ============================================================
function update_hL_from_alpha_frac_fast!(
    hL_vec::Vector{Float64},
    alpha_frac_vec::Vector{Float64},
    D::Float64,
    table::FastAlphaFracTable
)
    @assert length(hL_vec) == length(alpha_frac_vec)

    @inbounds for i in eachindex(alpha_frac_vec)
        hL_vec[i] = compute_hL_from_alpha_frac_fast(alpha_frac_vec[i], D, table)
    end

    return nothing
end

# ============================================================
# Atualização vetorial - D variável ao longo da malha
# ============================================================
function update_hL_from_alpha_frac_fast!(
    hL_vec::Vector{Float64},
    alpha_frac_vec::Vector{Float64},
    D_vec::Vector{Float64},
    table::FastAlphaFracTable
)
    @assert length(hL_vec) == length(alpha_frac_vec) == length(D_vec)

    @inbounds for i in eachindex(alpha_frac_vec)
        hL_vec[i] = compute_hL_from_alpha_frac_fast(alpha_frac_vec[i], D_vec[i], table)
    end

    return nothing
end