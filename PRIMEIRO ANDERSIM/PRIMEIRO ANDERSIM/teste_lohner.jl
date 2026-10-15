
# function lohner_indicator(U; epsilon=1e-8)
#     N = length(U)
#     E = zeros(N)

#     for i in 2:(N-1)
#         numerator = abs(U[i+1] - 2*U[i] + U[i-1])

#         denominator =
#             abs(U[i+1] - U[i]) +
#             abs(U[i] - U[i-1]) +
#             epsilon * (
#                 abs(U[i+1]) +
#                 2*abs(U[i]) +
#                 abs(U[i-1])
#             )

#         E[i] = numerator / max(denominator, eps(Float64))
#     end

#     return E
# end

# U = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

# E = lohner_indicator(U)

# println("U = ", U)
# println("Indicador de Löhner = ", E)
# println("Maior indicador = ", maximum(E))



# ============================================================
# TESTE DO INDICADOR DE LÖHNER
# Objetivo:
# Calcular e visualizar um indicador de variação local
# para identificar regiões candidatas ao refinamento da malha.
# ============================================================

using Plots

# ------------------------------------------------------------
# FUNÇÃO: lohner_indicator
#
# Entrada:
#   U       -> valores de uma variável nas células da malha.
#   epsilon -> pequeno parâmetro para evitar instabilidades
#              no denominador quando os valores são pequenos.
#
# Saída:
#   E       -> indicador de Löhner em cada célula.
#
# Interpretação:
#   Valores maiores indicam regiões onde a solução apresenta
#   mudanças locais acentuadas, que podem exigir refinamento.
# ------------------------------------------------------------
function lohner_indicator(U; epsilon=1e-8)

    # Número total de células representadas no vetor U.
    N = length(U)

    # Inicializa o indicador com zero em todas as células.
    E = zeros(N)

    # Calcula o indicador nas células internas.
    # As extremidades são ignoradas porque não possuem
    # vizinhos dos dois lados.
    for i in 2:(N-1)

        # Segunda diferença: mede a mudança da inclinação
        # da solução entre três células consecutivas.
        numerator = abs(U[i+1] - 2*U[i] + U[i-1])

        # Denominador que normaliza a segunda diferença
        # pelas variações entre células vizinhas.
        denominator =
            abs(U[i+1] - U[i]) +
            abs(U[i] - U[i-1]) +
            epsilon * (
                abs(U[i+1]) +
                2*abs(U[i]) +
                abs(U[i-1])
            )

        # Calcula o indicador, protegendo contra divisão
        # por um denominador extremamente pequeno.
        E[i] = numerator / max(denominator, eps(Float64))
    end

    # Devolve o vetor com o indicador calculado.
    return E
end

# ------------------------------------------------------------
# TESTE COM UMA SOLUÇÃO SIMPLES
#
# A solução muda de 0 para 1.
# Essa transição representa uma descontinuidade idealizada.
# ------------------------------------------------------------

# Valores da variável em seis células consecutivas.
U = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

# Calcula o indicador de Löhner para a solução de teste.
E = lohner_indicator(U)

# Imprime os resultados no terminal.
println("Valores de U = ", U)
println("Indicador de Löhner = ", E)
println("Maior indicador = ", maximum(E))

# ------------------------------------------------------------
# VISUALIZAÇÃO
#
# O primeiro gráfico mostra a solução U.
# O segundo mostra o indicador E.
# ------------------------------------------------------------

# Eixo horizontal: índice de cada célula.
celulas = 1:length(U)

# Gráfico da solução usada no teste.
grafico_solucao = plot(
    celulas, U;
    marker=:circle,
    linewidth=2,
    xlabel="Índice da célula",
    ylabel="U",
    title="Solução de teste",
    label="U",
    legend=:topleft
)

# Gráfico do indicador de Löhner.
grafico_indicador = plot(
    celulas, E;
    marker=:circle,
    linewidth=2,
    xlabel="Índice da célula",
    ylabel="Indicador de Löhner",
    title="Indicador de Löhner",
    label="E",
    legend=:topleft
)

# Combina os dois gráficos em uma única figura.
figura = plot(
    grafico_solucao,
    grafico_indicador;
    layout=(2, 1),
    size=(800, 600)
)

# Exibe a figura na janela de gráficos.
display(figura)

# Salva a figura na pasta atual do projeto.
savefig(figura, "teste_lohner.png")

println("Gráfico salvo em teste_lohner.png")
