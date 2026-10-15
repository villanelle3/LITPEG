############################## AnderSim v.1 ####################################
# Escoamento bifásico (gás-líquido) 1D transiente, incompressível, isotérmico,
# sem efeitos difusivos. Escrito inteiramente em modo serial (sem paralelismo).
############################ By Anderson Viana #################################
#
# ============================================================================
# OBSERVAÇÃO GERAL SOBRE A NOMENCLATURA USADA NESTE CÓDIGO
# ============================================================================
# O código usa três "camadas" de variáveis para cada campo (uL, uG, αL, αG, P):
#
#   - SEM sufixo  (ex.: uL)  → valor CORRIGIDO na iteração atual
#   - Sufixo 0    (ex.: uL0) → valor no INSTANTE DE TEMPO anterior (t - Δt)
#   - Sufixo Ø    (ex.: uLØ) → valor "chutado"/estimado na iteração anterior
#                              (usado como linearização dentro do laço interno)
#
# Além disso há as variáveis de correção (uLδ, uGδ, Pδ), que armazenam
# apenas o INCREMENTO calculado pela equação de pressão.
#
# Isso é a base de um esquema tipo SIMPLE (Semi-Implicit Method for
# Pressure-Linked Equations), muito comum em CFD de escoamentos multifásicos.
# ============================================================================

# ----------------------------------------------------------------------------
# Inclusão do solver TDMA (Thomas Algorithm) — resolve sistemas tridiagonais
# ----------------------------------------------------------------------------
include("./TDMA.jl")

# ----------------------------------------------------------------------------
# Bibliotecas utilizadas
# ----------------------------------------------------------------------------
using Printf          # Formatação de saída (printf)
using Plots           # Geração de gráficos
using CSV             # Escrita de arquivos CSV
using DataFrames      # Estrutura de dados tabular (DataFrame)
using LinearAlgebra   # Função norm() para cálculo de erros

# ============================================================================
# 1) LEITURA DO ARQUIVO DE ENTRADA (Inputs.txt)
# ============================================================================
n_dados = 16                       # Número de parâmetros de entrada
dados = zeros(n_dados)             # Vetor que armazenará os 16 valores
texto = readlines("Inputs.txt")    # Lê todas as linhas do arquivo

# O arquivo tem formato: linha ímpar = nome do parâmetro, linha par = valor.
# Portanto, percorremos apenas as linhas pares (i = 2, 4, 6, ...).
for i = 2:2:(2*n_dados)
    aux = texto[i]                 # Linha contendo o valor numérico (string)
    aux2 = Int(i / 2)              # Índice correspondente no vetor 'dados'
    global dados[aux2] = parse(Float64, aux)  # Converte string → Float64
end

# ============================================================================
# 2) VARIÁVEIS DE ENTRADA (mapeamento do vetor 'dados' para nomes)
# ============================================================================
Comp        = dados[1]    # Comprimento total do domínio (m)
D           = dados[2]    # Diâmetro interno do tubo (m)
θ           = dados[3]    # Ângulo de inclinação em relação à horizontal (rad)
n_segmentos = dados[4]    # Número de segmentos (malha uniforme)
tempo       = dados[5]    # Tempo total de simulação (s)
uL_inlet    = dados[6]    # Velocidade do líquido na entrada (m/s)
uG_inlet    = dados[7]    # Velocidade do gás na entrada (m/s)
ρL          = dados[8]    # Massa específica do líquido (kg/m³)
ρG          = dados[9]    # Massa específica do gás (kg/m³)
P_out       = dados[10]   # Pressão na saída do domínio (Pa)
αG_inlet    = dados[11]   # Fração volumétrica de gás na entrada (-)
αG_i        = dados[12]   # Fração volumétrica de gás inicial no domínio (-)
CFL         = dados[13]   # Número de Courant-Friedrichs-Lewy (estabilidade)
Tol_P       = dados[14]   # Tolerância de convergência para pressão
Tol_u       = dados[15]   # Tolerância de convergência para velocidades
Tol_α       = dados[16]   # Tolerância de convergência para frações

# ============================================================================
# 3) PARÂMETROS DERIVADOS E GEOMETRIA DA MALHA
# ============================================================================
N = Int(n_segmentos)                      # Número de nós da malha principal
A = (pi * D^2) / 4                        # Área da seção transversal (m²) — constante
αL_inlet = 1 - αG_inlet                   # Fração de líquido na entrada
αL_i     = 1 - αG_i                       # Fração de líquido inicial
u_inlet  = [uL_inlet, uG_inlet]           # Vetor com velocidades de entrada
Δx = Comp / n_segmentos                   # Tamanho do passo espacial (m)
aux = maximum(u_inlet)                    # Maior velocidade de entrada (para Δt inicial)
Δt = (CFL * Δx) / aux                     # Passo de tempo inicial (CFL)

# ============================================================================
# 4) ALOCAÇÃO DAS VARIÁVEIS PRINCIPAIS
# ============================================================================
# Campos no instante ATUAL
P   = zeros(N)        # Pressão nos N nós
αL  = zeros(N)        # Fração volumétrica de líquido
αG  = zeros(N)        # Fração volumétrica de gás
uL  = zeros(N+1)      # Velocidade do líquido nas N+1 faces
uG  = zeros(N+1)      # Velocidade do gás nas N+1 faces

# Campos no instante ANTERIOR (t - Δt) — usados nos termos transientes
P0   = zeros(N)
αL0  = zeros(N)
αG0  = zeros(N)
uL0  = zeros(N+1)
uG0  = zeros(N+1)

# Campos "chutados" na iteração anterior do laço interno (SIMPLE)
PØ   = zeros(N)
αLØ  = zeros(N)
αGØ  = zeros(N)
uLØ  = zeros(N+1)
uGØ  = zeros(N+1)

# ============================================================================
# 5) CONDIÇÕES INICIAIS
# ============================================================================
uL[:] .= uL_inlet        # Velocidade do líquido = valor de entrada
uG[:] .= uG_inlet        # Velocidade do gás = valor de entrada
αL[:] .= αL_i            # Fração de líquido uniforme = valor inicial
αG[:] .= αG_i            # Fração de gás uniforme = valor inicial
P[:]  .= P_out           # Pressão uniforme = pressão de saída

# ============================================================================
# 6) MAPA DE VALORES (inicialização das três camadas de variáveis)
# ============================================================================
# Camada "tempo anterior"
uL0[:] .= uL_inlet
uG0[:] .= uG_inlet
αL0[:] .= αL_i
αG0[:] .= αG_i
P0[:]  .= P_out

# Camada "chute" (iteração anterior do laço interno)
uLØ[:] .= uL_inlet
uGØ[:] .= uG_inlet
αLØ[:] .= αL_i
αGØ[:] .= αG_i
PØ[:]  .= P_out

# ============================================================================
# 7) VARIÁVEIS DE CORREÇÃO (incrementos)
# ============================================================================
uLδ = zeros(N+1)   # Correção da velocidade do líquido
uGδ = zeros(N+1)   # Correção da velocidade do gás
Pδ  = zeros(N)     # Correção da pressão (resolvida pela eq. de Poisson)

# ============================================================================
# 8) SISTEMAS LINEARES — VELOCIDADES (matrizes tridiagonais)
# ============================================================================
# Convenção: d = diagonal principal, du = diagonal superior, dl = diagonal inferior
#            b = lado direito (RHS)
d_L  = zeros(N+1)  # Diagonal principal    — líquido
du_L = zeros(N+1)  # Diagonal superior     — líquido
dl_L = zeros(N+1)  # Diagonal inferior     — líquido
b_L  = zeros(N+1)  # RHS                   — líquido

d_G  = zeros(N+1)  # Diagonal principal    — gás
du_G = zeros(N+1)  # Diagonal superior     — gás
dl_G = zeros(N+1)  # Diagonal inferior     — gás
b_G  = zeros(N+1)  # RHS                   — gás

# ----------------------------------------------------------------------------
# Condições de contorno para as VELOCIDADES
# ----------------------------------------------------------------------------
# Face 1 (entrada): Dirichlet — velocidade prescrita
d_L[1]  = 1.0
du_L[1] = 0.0
dl_L[1] = 0.0
b_L[1]  = uL_inlet
d_G[1]  = 1.0
du_G[1] = 0.0
dl_G[1] = 0.0
b_G[1]  = uG_inlet

# Face N+1 (saída): Neumann — gradiente nulo (u[N+1] - u[N] = 0)
d_L[N+1]  = 1.0
du_L[N+1] = 0.0
dl_L[N+1] = -1.0
b_L[N+1]  = 0
d_G[N+1]  = 1.0
du_G[N+1] = 0.0
dl_G[N+1] = -1.0
b_G[N+1]  = 0

# ============================================================================
# 9) SISTEMAS LINEARES — FRAÇÕES VOLUMÉTRICAS (αL, αG)
# ============================================================================
d_Lα  = zeros(N)
du_Lα = zeros(N)
dl_Lα = zeros(N)
b_Lα  = zeros(N)

d_Gα  = zeros(N)
du_Gα = zeros(N)
dl_Gα = zeros(N)
b_Gα  = zeros(N)

# ----------------------------------------------------------------------------
# Condições de contorno para as FRAÇÕES
# ----------------------------------------------------------------------------
# Nó 1 (entrada): Dirichlet — fração prescrita
d_Lα[1]  = 1.0
du_Lα[1] = 0.0
dl_Lα[1] = 0.0
b_Lα[1]  = αL_inlet
d_Gα[1]  = 1.0
du_Gα[1] = 0.0
dl_Gα[1] = 0.0
b_Gα[1]  = αG_inlet

# Nó N (saída): Neumann — gradiente nulo (α[N] - α[N-1] = 0)
d_Lα[N]  = 1.0
du_Lα[N] = 0.0
dl_Lα[N] = -1.0
b_Lα[N]  = 0
d_Gα[N]  = 1.0
du_Gα[N] = 0.0
dl_Gα[N] = -1.0
b_Gα[N]  = 0.0

# ============================================================================
# 10) SISTEMA LINEAR — PRESSÃO (equação de Poisson para correção)
# ============================================================================
d  = zeros(N)   # Diagonal principal
du = zeros(N)   # Diagonal superior
dl = zeros(N)   # Diagonal inferior
b  = zeros(N)   # RHS

# ----------------------------------------------------------------------------
# Condições de contorno para a PRESSÃO (a variável resolvida é Pδ, a correção)
# ----------------------------------------------------------------------------
# Nó 1 (entrada): gradiente nulo → Pδ[1] - Pδ[2] = 0
d[1]  = 1.0
du[1] = -1.0
dl[1] = 0.0
b[1]  = 0.0

# Nó N (saída): Dirichlet → Pδ[N] = 0 (a pressão de saída já é conhecida)
d[N]  = 1.0
du[N] = 0.0
dl[N] = 0.0
b[N]  = 0.0

# ============================================================================
# 11) LAÇO PRINCIPAL NO TEMPO
# ============================================================================
t = 0
while t < tempo

    # ------------------------------------------------------------------------
    # 11.1) Cálculo do passo de tempo Δt via condição CFL
    # ------------------------------------------------------------------------
    aux1 = maximum(abs.(uL[:]))    # Maior |uL| atual
    aux2 = maximum(abs.(uG[:]))    # Maior |uG| atual
    #aux1 += cL                     # (comentado) velocidade do som — não usado
    #aux2 += cG                     # (comentado) velocidade do som — não usado
    u = [aux1, aux2]
    wavevel = maximum(u)           # Velocidade de onda máxima no domínio
    global Δt = (CFL * Δx) / (wavevel)
    global t
    global Δt = minimum([Δt, (tempo - t)])   # Não ultrapassar o tempo final
    global t = t + Δt                        # Avança o tempo
    print(t)
    print("\n")

    # ------------------------------------------------------------------------
    # 11.2) Declaração 'global' (necessária em scripts Julia com while)
    # ------------------------------------------------------------------------
    global uL0;  global uG0;  global αL0;  global αG0;  global P0
    global uLØ;  global uGØ;  global αLØ;  global αGØ;  global PØ
    global uL;   global uG;   global αL;   global αG;   global P

    # ------------------------------------------------------------------------
    # 11.3) Correlações empíricas (ainda NÃO implementadas!)
    # ------------------------------------------------------------------------
    τ_L = 0.0      # Tensão de cisalhamento parede-líquido     — a implementar
    τ_G = 0.0      # Tensão de cisalhamento parede-gás        — a implementar
    Γ_G = 0.0      # Taxa de transferência de massa interfacial (gás) — a implementar
    Γ_L = -Γ_G     # Conservação de massa: Γ_L = -Γ_G

    cont = 0       # Contador (não utilizado)

    # ------------------------------------------------------------------------
    # 11.4) Inicialização dos erros (força entrar no laço interno)
    # ------------------------------------------------------------------------
    e_uL = 1
    e_uG = 1
    e_αL = 1
    e_αG = 1
    e_P  = 1

    aux_uLØ = zeros(N+1)   # Armazena uLØ antes da atualização (para erros)
    aux_uGØ = zeros(N+1)   # Armazena uGØ antes da atualização (para erros)

    cont1 = 0

    # ========================================================================
    # 12) LAÇO INTERNO DE ACOPLAMENTO (tipo SIMPLE / Picard)
    # ========================================================================
    while (e_uL>Tol_u || e_uG>Tol_u || e_αL>Tol_α || e_αG>Tol_α || e_P>Tol_P) && cont1 < 1000

        cont1 += 1

        # --------------------------------------------------------------------
        # 12.1) EQUAÇÃO DE QUANTIDADE DE MOVIMENTO — VELOCIDADES (uL, uG)
        # --------------------------------------------------------------------

        # Reaplica as condições de contorno de saída (Neumann — gradiente nulo)
        d_L[N+1]  = 1.0
        du_L[N+1] = 0.0
        dl_L[N+1] = -1.0
        b_L[N+1]  = 0
        d_G[N+1]  = 1.0
        du_G[N+1] = 0.0
        dl_G[N+1] = -1.0
        b_G[N+1]  = 0

        # Percorre as faces internas (i = 1 .. N-1)
        for i=1:(N-1)

            # ----------------------------------------------------------------
            # Velocidade relativa interfacial (para o termo de arraste Γ·ui)
            # Ψ = 0 → usa uLØ (se Γ_G >= 0);  Ψ = 1 → usa uGØ (se Γ_G < 0)
            # ----------------------------------------------------------------
            if Γ_G >= 0
                Ψ = 0
            else
                Ψ = 1
            end
            ui = Ψ * uGØ[i+1] + (1 - Ψ) * uLØ[i+1]

            # ----------------------------------------------------------------
            # Coeficientes pós-agrupamento (discretização upwind)
            # ----------------------------------------------------------------
            # Termos transientes (massa × 1/Δt) — média entre nós i e i+1
            mt_L = (ρL * αL0[i] + ρL * αL0[i+1]) / (2 * Δt)
            mt_G = (ρG * αG0[i] + ρG * αG0[i+1]) / (2 * Δt)

            # Termos de gradiente de pressão (α/Δx) — média entre nós
            mx_L = (αLØ[i] + αLØ[i+1]) / (2 * Δx)
            mx_G = (αGØ[i] + αGØ[i+1]) / (2 * Δx)

            # Termos gravitacionais (ρ·α·g·sinθ) — média entre nós
            mg_L = 9.81*(ρL*αLØ[i] + ρL*αLØ[i+1])*sin(θ)/2
            mg_G = 9.81*(ρG*αGØ[i] + ρG*αGØ[i+1])*sin(θ)/2

            # Fluxos de massa nas faces leste (e) e oeste (w) — líquido
            G_Les = ρL * αLØ[i+1] * uLØ[i+1]   # face leste
            G_Lws = ρL * αLØ[i]   * uLØ[i]     # face oeste

            # Fluxos de massa — gás (upwind dependendo do sinal de uG)
            if uGØ[i+1] >= 0
                G_Ges = ρG * αGØ[i+1] * uGØ[i+1]
                G_Gws = ρG * αGØ[i]   * uGØ[i]
            else
                G_Ges = ρG * αGØ[i+1] * uGØ[i+2]   # upwind "a jusante"
                G_Gws = ρG * αGØ[i]   * uGØ[i+1]
            end

            # ----------------------------------------------------------------
            # Montagem do sistema linear para o LÍQUIDO
            #   d_L·uL[i+1] + dl_L·uL[i] = b_L
            # ----------------------------------------------------------------
            d_L[i+1]  = mt_L + G_Lws / Δx + Γ_L
            du_L[i+1] = 0.0
            dl_L[i+1] = -G_Lws / Δx
            b_L[i+1]  = mt_L*uL0[i+1] - mx_L*(PØ[i+1]-PØ[i]) - mg_L - Γ_L + Γ_L*ui

            # ----------------------------------------------------------------
            # Montagem do sistema linear para o GÁS (upwind conforme sinal)
            # ----------------------------------------------------------------
            if uGØ[i+1] >= 0
                d_G[i+1]  = mt_G + (G_Gws/Δx) + Γ_G
                du_G[i+1] = 0.0
                dl_G[i+1] = -G_Gws/Δx
            else
                d_G[i+1]  = mt_G + Γ_G - G_Ges/Δx
                du_G[i+1] = G_Ges/Δx
                dl_G[i+1] = 0.0
            end
            b_G[i+1] = mt_G*uG0[i+1] - mx_G*(PØ[i+1]-PØ[i]) - mg_G - Γ_G + Γ_G*ui
        end

        # Salva os valores antigos para o cálculo dos erros
        aux_uLØ[:] = uLØ[:]
        aux_uGØ[:] = uGØ[:]

        # Resolve os dois sistemas tridiagonais (TDMA)
        uLØ = TDMA(d_L, du_L, dl_L, b_L)
        uGØ = TDMA(d_G, du_G, dl_G, b_G)

        # Extrapolação linear para a face de saída (N+1)
        uLØ[N+1] = 2*uLØ[N] - uLØ[N-1]
        uGØ[N+1] = 2*uGØ[N] - uGØ[N-1]

        # --------------------------------------------------------------------
        # 12.2) EQUAÇÃO DE PRESSÃO (Poisson para correção Pδ)
        # --------------------------------------------------------------------

        # Condição de contorno na saída (Dirichlet — Pδ = 0)
        d[N]  = 1.0
        du[N] = 0.0
        dl[N] = 0.0
        b[N]  = 0.0

        # Percorre os nós internos (i = 2 .. N-1)
        for i = 2:(N-1)

            # Transmissibilidades (coeficientes tipo difusão) nas faces e/w
            d_Le = (Δt / Δx) * (αLØ[i] + αLØ[i+1]) / (ρL * αLØ[i] + ρL * αLØ[i+1])
            d_Lw = (Δt / Δx) * (αLØ[i] + αLØ[i-1]) / (ρL * αLØ[i] + ρL * αLØ[i-1])
            d_Ge = (Δt / Δx) * (αGØ[i] + αGØ[i+1]) / (ρG * αGØ[i] + ρG * αGØ[i+1])
            d_Gw = (Δt / Δx) * (αGØ[i] + αGØ[i-1]) / (ρG * αGØ[i] + ρG * αGØ[i-1])

            # Fluxos volumétricos (α·u) nas faces — líquido
            G_Le = αLØ[i]   * uLØ[i+1]
            G_Lw = αLØ[i-1] * uLØ[i]

            # Fluxos volumétricos — gás (upwind conforme sinal de uG)
            if uG[i] >= 0
                G_Ge    = αGØ[i]   * uGØ[i+1]
                G_Gw    = αGØ[i-1] * uGØ[i]
                alfaG_e = αGØ[i]
                alfaG_w = αGØ[i-1]
            else
                G_Ge    = αGØ[i+1] * uGØ[i+1]
                G_Gw    = αGØ[i]   * uGØ[i]
                alfaG_e = αGØ[i+1]
                alfaG_w = αGØ[i]
            end

            # Montagem da equação de Poisson para a correção de pressão
            d[i]  = (alfaG_e * d_Ge + alfaG_w * d_Gw + αLØ[i] * d_Le +
                     αLØ[i-1] * d_Lw)
            du[i] = -(alfaG_e * d_Ge + αLØ[i] * d_Le)
            dl[i] = -alfaG_w * d_Gw - αLØ[i-1] * d_Lw

            # RHS: balanço de massa (transiente + advecção) → resíduo a zerar
            b[i] = (-(Δx * αGØ[i]) / (Δt) + (Δx * αG0[i]) / (Δt) -
                    (Δx * αLØ[i]) / (Δt) + (Δx * αL0[i]) / (Δt)
                    - G_Ge + G_Gw - G_Le + G_Lw)
        end

        # Resolve o sistema tridiagonal da pressão
        Pδ = TDMA(d, du, dl, b)

        # Extrapolação para o nó 1 (entrada) — gradiente nulo
        Pδ[1] = 2*Pδ[2] - Pδ[3]

        # --------------------------------------------------------------------
        # 12.3) CORREÇÕES — aplica Pδ para corrigir velocidades
        # --------------------------------------------------------------------
        for i = 1:(N-1)
            de_L = (Δt/Δx)*(αLØ[i] + αLØ[i+1])/(ρL*αLØ[i] + ρL*αLØ[i+1])
            de_G = (Δt/Δx)*(αGØ[i] + αGØ[i+1])/(ρG*αGØ[i] + ρG*αGØ[i+1])
            uLδ[i+1] = de_L * (Pδ[i] - Pδ[i+1])
            uGδ[i+1] = de_G * (Pδ[i] - Pδ[i+1])
        end

        # Atualiza pressão e velocidades corrigidas
        P[:]  = PØ[:] + Pδ[:]
        uL    = uLØ + uLδ
        uG    = uGØ + uGδ

        # Reaplica condição de contorno na face de saída
        uL[N+1] = uL[N]
        uG[N+1] = uG[N]

        # --------------------------------------------------------------------
        # 12.4) EQUAÇÕES DE CONSERVAÇÃO DE MASSA — FRAÇÕES (αL, αG)
        # --------------------------------------------------------------------

        # Condição de contorno de saída (Neumann — gradiente nulo)
        d_Lα[N]  = 1.0
        du_Lα[N] = 0.0
        dl_Lα[N] = -1.0
        b_Lα[N]  = 0
        d_Gα[N]  = 1.0
        du_Gα[N] = 0.0
        dl_Gα[N] = -1.0
        b_Gα[N]  = 0.0

        # Percorre os nós internos (i = 2 .. N-1)
        for i = 2:(N-1)

            # ----------------------------------------------------------------
            # Líquido: transiente + advecção upwind
            # ----------------------------------------------------------------
            d_Lα[i]  = ρL / Δt + ρL * uL[i+1] / Δx
            du_Lα[i] = 0.0
            dl_Lα[i] = -ρL * uL[i] / Δx
            b_Lα[i]  = ρL * αL0[i] / Δt + Γ_L

            # ----------------------------------------------------------------
            # Gás: upwind conforme sinal de uG
            # ----------------------------------------------------------------
            if uG[i] >= 0
                d_Gα[i]  = ρG / Δt + ρG * uG[i+1] / Δx
                du_Gα[i] = 0.0
                dl_Gα[i] = -ρG * uG[i] / Δx
                b_Gα[i]  = ρG * αG0[i] / Δt + Γ_G
            else
                d_Gα[i]  = ρG / Δt - ρG * uG[i] / Δx
                du_Gα[i] = ρG * uG[i+1] / Δx
                dl_Gα[i] = 0.0
                b_Gα[i]  = ρG * αG0[i] / Δt + Γ_G
            end
        end

        # Resolve os dois sistemas (frações "brutas")
        αLØ = TDMA(d_Lα, du_Lα, dl_Lα, b_Lα)
        αGØ = TDMA(d_Gα, du_Gα, dl_Gα, b_Gα)

        # --------------------------------------------------------------------
        # 12.5) NORMALIZAÇÃO DAS FRAÇÕES (garante αL + αG = 1)
        # --------------------------------------------------------------------
        for j = 2:N
            αL[j] = αLØ[j] / (αLØ[j] + αGØ[j])
            αG[j] = 1.0 - αL[j]
        end

        # --------------------------------------------------------------------
        # 12.6) CÁLCULO DOS ERROS (norma relativa L2)
        # --------------------------------------------------------------------
        e_uL = norm(uL-uLØ)/norm(uLØ)
        e_uG = norm(uG-uGØ)/norm(uGØ)
        e_αL = norm(αL-αLØ)/norm(αLØ)
        e_αG = norm(αG-αGØ)/norm(αGØ)
        e_P  = norm(P-PØ)/norm(PØ)

        # --------------------------------------------------------------------
        # 12.7) ATUALIZAÇÃO DOS VALORES "CHUTADOS" PARA A PRÓXIMA ITERAÇÃO
        # --------------------------------------------------------------------
        uLØ[:] = uL[:]
        uGØ[:] = uG[:]
        αLØ[:] = αL[:]
        αGØ[:] = αG[:]
        PØ[:]  = P[:]
    end

    # ------------------------------------------------------------------------
    # 13) ATUALIZAÇÃO TEMPORAL — o instante atual vira o "anterior"
    # ------------------------------------------------------------------------
    uL0[:] = uL[:]
    uG0[:] = uG[:]
    αL0[:] = αL[:]
    αG0[:] = αG[:]
    P0[:]  = P[:]
end

# ============================================================================
# 14) SAÍDA DE DADOS
# ============================================================================
# Cria um DataFrame com os resultados finais e exporta para CSV
df = DataFrame(Pr=P, UL=uL[1:N], UG=uG[1:N], AlfaL=αL, AlfaG=αG)
CSV.write("export.csv", df)

# ============================================================================
# 15) GRÁFICO FINAL — perfil de fração de gás ao longo do domínio
# ============================================================================
aux2 = Comp/(N-1)              # Passo espacial usado para o eixo x do gráfico
plot(0:aux2:Comp, αG)          # Plota αG em função da posição