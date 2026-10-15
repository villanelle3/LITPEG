using Glob, DelimitedFiles, LaTeXStrings

"""
    quadPointsAndWeights(elementType = 1, p = 2, quadType = "legendre"; T = Float64)

Reads the integration point location and the Gauss-quadrature weight according to the (integer) polynomial order `p`,
and the (string) `quadratureType`. The only supported `elementType` is 1.
The `quadratureType` takes can take the value:
- "legendre"
- "equidistant"
- "lobatto"
Note that the Gauss-Lobatto flux reconstruction is still not implemented.
"""
function quadPointsAndWeights(elementType = 1, p = 2, quadType = "lobatto"; T = Float64)
    projectPath(parts...) = normpath(joinpath(@__DIR__, parts...))
    src = projectPath("../quadratures/")
    for (root, dirs, _) in walkdir(src)
        for dir in dirs
            if occursin(string(elementType), dir)
                file = glob("*" * quadType * "*" * string(p) * "*", joinpath(root, dir))[1]
                data = readdlm(file, '\t', T, skipstart = 1)
                return data[:, 1] , data[:, 2]
            end
        end
    end
end