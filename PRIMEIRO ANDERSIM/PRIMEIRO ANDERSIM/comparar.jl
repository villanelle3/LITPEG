# using CSV, DataFrames, Plots

# ref = CSV.read("ref_original.csv", DataFrame)
# amr = CSV.read("malha_fixa_L1.csv", DataFrame)

# dx = 12.0 / nrow(ref)
# xref = [(i-0.5)*dx for i in 1:nrow(ref)]

# p = plot(xref, ref.AlfaG,
#     label="uniforme N=256",
#     marker=:circle, ms=2)

# plot!(p, amr.x, amr.AlfaG,
#     label="malha fixa refinada",
#     marker=:circle, ms=2)

# xlims!(p, 5, 8)

# savefig(p, "comparacao.png")
# println("Gráfico salvo em comparacao.png")

# using CSV, DataFrames, Plots

# ref  = CSV.read("ref_original.csv", DataFrame)
# amr  = CSV.read("malha_fixa_0a7.csv", DataFrame)
# fino = CSV.read("uniforme_N512.csv", DataFrame)

# fino.xr = round.(fino.x, digits=8)
# amr.xr  = round.(amr.x, digits=8)

# m = innerjoin(fino, amr, on=:xr, makeunique=true)
# zona = (m.xr .>= 5.6) .& (m.xr .<= 7.4)

# println("pontos em comum: ", nrow(m))
# println("max |dAlfaG| na zona refinada: ",
#         maximum(abs.(m.AlfaG[zona] .- m.AlfaG_1[zona])))
# println("max |dAlfaG| no duto todo:     ",
#         maximum(abs.(m.AlfaG .- m.AlfaG_1)))

# dx = 12.0 / nrow(ref)

# p = plot(fino.x, fino.AlfaG,
#     label="uniforme 512")

# plot!(p, [(i-0.5)*dx for i in 1:nrow(ref)], ref.AlfaG,
#     label="uniforme 256")

# plot!(p, amr.x, amr.AlfaG,
#     label="refinada fixa", ls=:dash)

# xlims!(p, 5, 8)

# # Salvar a figura
# savefig(p, "comparacao_malhas.png")
# savefig(p, "comparacao_malhas.pdf")

using CSV, DataFrames, Plots

ref  = CSV.read("ref_original.csv", DataFrame)
amr  = CSV.read("malha_fixa_0a7.csv", DataFrame)
fino = CSV.read("uniforme_N512.csv", DataFrame)

fino.xr = round.(fino.x, digits=8)
amr.xr  = round.(amr.x, digits=8)

m = innerjoin(fino, amr, on=:xr, makeunique=true)

# Interpolação da solução uniforme N=256
function interp(xs, ys, x)
    i = clamp(searchsortedlast(xs, x), 1, length(xs)-1)
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return ys[i]*(1-t) + ys[i+1]*t
end

dx = 12.0 / nrow(ref)
x256 = [(i-0.5)*dx for i in 1:nrow(ref)]

m.a256 = [interp(x256, ref.AlfaG, x) for x in m.xr]

# Comparação no intervalo 0,1 <= x <= 7,0
sel = (m.xr .>= 0.1) .& (m.xr .<= 7.0)

println("pontos em comum: ", nrow(m))
println("max |refinada - 512| (0,1 <= x <= 7): ",
        maximum(abs.(m.AlfaG_1[sel] .- m.AlfaG[sel])))

println("max |256 - 512| (0,1 <= x <= 7): ",
        maximum(abs.(m.a256[sel] .- m.AlfaG[sel])))

# Gráfico comparativo
p = plot(fino.x, fino.AlfaG,
    label="uniforme 512",
    xlabel="x (m)",
    ylabel="AlfaG")

plot!(p, x256, ref.AlfaG,
    label="uniforme 256")

plot!(p, amr.x, amr.AlfaG,
    label="refinada 0 a 7,5 m",
    ls=:dash)

xlims!(p, 5, 8)

# Salvar as figuras
savefig(p, "comparacao_malhas_0a7.png")
savefig(p, "comparacao_malhas_0a7.pdf")