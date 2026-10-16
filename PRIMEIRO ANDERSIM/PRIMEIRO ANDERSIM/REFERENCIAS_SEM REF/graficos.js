
const fs = require("fs");
const { parse } = require("csv-parse/sync");
const { ChartJSNodeCanvas } = require("chartjs-node-canvas");

async function main() {
    // 1. Ler o arquivo CSV
    const csv = fs.readFileSync("referencia_256.csv", "utf8");

    const dados = parse(csv, {
        columns: true,
        skip_empty_lines: true,
        trim: true
    });

    // 2. Converter os valores para números
    const pontos = dados.map((d, i) => ({
        x: i * 12 / (dados.length - 1),
        AlfaG: Number(d.AlfaG),
        UL: Number(d.UL)
    }));

    // 3. Obter os valores de entrada e saída
    const entrada = pontos[0];
    const saida = pontos[pontos.length - 1];

    // 4. Identificar a maior variação de AlfaG
    let indiceFrente = 0;
    let maiorGradiente = 0;

    for (let i = 0; i < pontos.length - 1; i++) {
        const dx = pontos[i + 1].x - pontos[i].x;

        const gradiente = Math.abs(
            (pontos[i + 1].AlfaG - pontos[i].AlfaG) / dx
        );

        if (gradiente > maiorGradiente) {
            maiorGradiente = gradiente;
            indiceFrente = i;
        }
    }

    // 5. Criar a ferramenta para gerar imagens
    const grafico = new ChartJSNodeCanvas({
        width: 1000,
        height: 500,
        backgroundColour: "white"
    });

    // 6. Função para gerar cada gráfico
    async function gerarGrafico(variavel, titulo, cor, arquivo) {
        const configuracao = {
            type: "line",

            data: {
                datasets: [{
                    label: variavel,
                    data: pontos.map(p => ({
                        x: p.x,
                        y: p[variavel]
                    })),
                    borderColor: cor,
                    borderWidth: 2,
                    pointRadius: 0,
                    tension: 0
                }]
            },

            options: {
                responsive: false,

                plugins: {
                    title: {
                        display: true,
                        text: titulo,
                        font: { size: 18 }
                    },
                    legend: {
                        display: true
                    }
                },

                scales: {
                    x: {
                        type: "linear",
                        min: 0,
                        max: 12,
                        title: {
                            display: true,
                            text: "Posição ao longo do duto (m)"
                        }
                    },

                    y: {
                        title: {
                            display: true,
                            text: variavel
                        }
                    }
                }
            }
        };

        const imagem = await grafico.renderToBuffer(configuracao);
        fs.writeFileSync(arquivo, imagem);

        console.log(`Gráfico salvo: ${arquivo}`);
    }

    // 7. Gerar os dois gráficos
    await gerarGrafico(
        "AlfaG",
        "Fração volumétrica de gás ao longo do duto",
        "blue",
        "AlfaG.png"
    );

    await gerarGrafico(
        "UL",
        "Velocidade do líquido ao longo do duto",
        "darkorange",
        "UL.png"
    );

    // 8. Imprimir os resultados
    console.log("\n========== RESULTADOS ==========");

    console.log(
        `AlfaG na entrada (x = 0 m): ${entrada.AlfaG.toFixed(6)}`
    );

    console.log(
        `AlfaG na saída (x = 12 m): ${saida.AlfaG.toFixed(6)}`
    );

    console.log(
        `Maior variação de AlfaG entre x = ${
            pontos[indiceFrente].x.toFixed(4)
        } m e x = ${
            pontos[indiceFrente + 1].x.toFixed(4)
        } m`
    );

    console.log(
        `Gradiente máximo de AlfaG: ${maiorGradiente.toFixed(6)} por metro`
    );

    console.log(
        `UL na entrada: ${entrada.UL.toFixed(6)} m/s`
    );

    console.log(
        `UL na saída: ${saida.UL.toFixed(6)} m/s`
    );

    console.log(
        `Entrada de UL igual a 10 m/s: ${
            Math.abs(entrada.UL - 10) < 1e-6 ? "SIM" : "NÃO"
        }`
    );
}

main().catch(erro => {
    console.error("Erro ao executar o programa:", erro.message);
});
