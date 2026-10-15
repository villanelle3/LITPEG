module AnderSim5EqWang
################################################################################
############################   AnderSim_v.Wang 5 Eq.  ##########################
############################## By Anderson Viana ###############################

using Printf
using Plots
using CSV
using DataFrames
using LinearAlgebra

include("functions/vandermonde_legendre.jl")
include("functions/TDMA.jl")
include("functions/Limiter.jl")
include("functions/hMLP.jl")
include("functions/minmod2.jl")
include("functions/hmlp2.jl")
include("functions/utils.jl")
include("functions/solver.jl")
include("functions/polynomials.jl")
include("functions/mesh.jl")
include("functions/flow.jl")
include("functions/frgrid.jl")
include("functions/FR.jl")
include("functions/MINMOD.jl")
include("functions/MUSCL.jl")
include("functions/closures_stratified.jl")
include("functions/hLCompute.jl")

export run_case

include("src/types_globals.jl")
include("src/io.jl")
include("src/allocation.jl")
include("src/initialization.jl")
include("src/faces.jl")
include("src/momentum.jl")
include("src/pressure.jl")
include("src/volume_fraction.jl")
include("src/energy.jl")
include("src/post.jl")
include("src/run_case.jl")

end # module
