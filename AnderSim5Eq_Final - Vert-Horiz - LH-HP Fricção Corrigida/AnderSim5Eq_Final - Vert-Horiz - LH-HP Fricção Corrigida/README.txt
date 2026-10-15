Estrutura:

- main.jl              -> módulo principal + includes
- run_main.jl          -> script simples para executar
- src/types_globals.jl -> structs leves e globais esperados pelos arquivos legados
- src/io.jl            -> leitura da entrada e construção da malha
- src/allocation.jl    -> alocação dos arrays
- src/initialization.jl-> inicialização global, do estado e do FR
- src/run_case.jl      -> laço principal da simulação

Como executar no Julia:

include("main.jl")
using .AnderSim5EqWang
run_case("Entrada.txt")

ou

include("run_main.jl")


Refinamento adicional:
- src/post.jl            -> exportação CSV + plot final
- src/faces.jl           -> reservado para extração futura
- src/momentum.jl        -> reservado para extração futura
- src/pressure.jl        -> reservado para extração futura
- src/volume_fraction.jl -> reservado para extração futura
- src/energy.jl          -> reservado para extração futura

