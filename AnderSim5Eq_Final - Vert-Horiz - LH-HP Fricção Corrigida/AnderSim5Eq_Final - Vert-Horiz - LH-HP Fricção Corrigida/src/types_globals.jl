Base.@kwdef struct InputData
    Comp::Float64
    D::Float64
    θ::Float64
    n_segmentos::Int
    tempo::Float64
    uL_inlet::Float64
    uG_inlet::Float64
    ρLi::Float64
    ρGi::Float64
    P_out::Float64
    αG_inlet::Float64
    αG_0::Float64
    T_intet::Float64
    T_0::Float64
    CFL::Float64
    Tol_P::Float64
    Tol_u::Float64
    Tol_α::Float64
    Cp_L::Float64
    Cp_G::Float64
    kf_L::Float64
    kf_G::Float64
    T_ext::Float64
    μL::Float64
    μG::Float64
    kf_pipe::Float64
    D_ext::Float64
    g::Float64
    cG::Float64
    cL::Float64
    βG::Float64
    βL::Float64
    γ::Float64
    h_ext::Float64
    metodo::Int
    Grau::Int
    C::Int
    energyflag::Int
    flowpattern::Int
    roughness::Float64
    estabmet::Int
end

Base.@kwdef struct GridData
    N::Int
    Δx::Float64
    LMIN::Float64
    LMAX::Float64
    A::Float64
    S::Float64
end

# ============================================================
# Globais esperados por flow.jl e afins
# ============================================================
global αG_inlet = 0.0
global αL_inlet = 0.0
global uL_inlet = 0.0
global uG_inlet = 0.0
global hL_intet = 0.0
global hLi = 0.0
global Grau = 0
global Λ = 1.0
global V = zeros(1,1)
global N = 0

# ============================================================
# Leitura da entrada
# ============================================================
