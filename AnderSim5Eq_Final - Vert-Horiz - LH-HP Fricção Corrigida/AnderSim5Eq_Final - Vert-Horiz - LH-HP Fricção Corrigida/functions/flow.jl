"""
    Flow{T <: AbstractFloat}

Derived type storing flow related quantities.

# Fields
- `uδD :: Array{T, 3}`. Discontinuous conserved quantities. Stores values at each [element, integration point, equation index].
- `wδD :: Array{T, 3}`. Discontinuous primitive quantities. Stores values at each [element, integration point, equation index].
- `fδD :: Array{T, 3}`. Discontinuous advective flux quantities. Stores values at [element, integration point, equation index].
- `qδD :: Array{T, 3}`. Discontinuous viscous flux quantities. Stores values at [element, integration point, equation index].
- `rhs :: Array{T, 3}`. RHS quantities (du/dt = RHS). Stores values at each [element, integration point, equation index].
- `uδI :: Matrix{T}`. Interaction (aka common, numerical) flux of conserved quantities. Stores values at each [face, equation index].
- `fδI :: Matrix{T}`. Interaction advective flux. Stores values at each [face, equation index].
- `qδI :: Matrix{T}`. Interaction viscous flux. Stores values at each [face, equation index].
- `uK :: Array{T, 3}`. Left/Right conserved quantities at the faces. Stores values at each [face, left/right side, equation index].
- `fK :: Array{T, 3}`. Left/Right advective flux at the faces. Stores values at each [face, left/right side, equation index].
- `qK :: Array{T, 3}`. Left/Right viscous flux at the faces. Stores values at each [face, left/right side, equation index].
- `wBC :: Array{T, 3}`. Primitive quantities stored at the physical boundaries (acting as BCs).
    Stores values at each [left/right boundary element, integration point, equation index].
- `δt :: Vector{T}`. Stores the time step for each update of the equations.
- `time :: Vector{T}`. Stores the simulation time for each update of the equations.
- `PDE :: Integer`. Flag defining the type of equations to be solved:
    1 ≡ Linear advection equation, 2 ≡ Burgers equation, 3 ≡ Viscous Burgers equation, 4 ≡ Euler equations.
- `BC :: Integer`. Flag defining the type of boundary conditions:
    0 ≡ Periodic, 1 ≡ Outlet.
- `γ :: AbstractFloat`. Ratio of specific heats.
- `ν :: AbstractFloat`. Viscosity.
# Arguments
- `N :: Integer`. Number of elements (excluding ghost elements).
- `P :: Integer`. Order of the elements.
- `PDE :: Integer`. (As above).
- `BC :: Integer`. (As above).
- `γ :: AbstractFloat`. (As above).
- `ν :: AbstractFloat`. (As above).
- `T :: Type`. Float precision type.
"""
struct Flow{T <: AbstractFloat}
    uδD :: Array{T, 3}
    wδD :: Array{T, 3}
    fδD :: Array{T, 3}
    qδD :: Array{T, 3}
    rhs :: Array{T, 3}
    uδI :: Matrix{T}
    fδI :: Matrix{T}
    qδI :: Matrix{T}
    uK :: Array{T, 3}
    fK :: Array{T, 3}
    qK :: Array{T, 3}
    wBC :: Array{T, 3}
    δt :: Vector{T}
    time :: Vector{T}
    PDE :: Integer
    BC :: Integer
    γ :: AbstractFloat
    ν :: AbstractFloat
    
"""
    Flow(N :: Integer, P :: Integer, PDE :: Integer, BC :: Integer;
        γ :: AbstractFloat = 1.4, ν :: AbstractFloat = 0.0, T = Float64)

Initializes the `Flow` struct.
"""
    function Flow(N :: Integer, P :: Integer, PDE :: Integer, BC :: Integer;
        γ :: AbstractFloat = 1.4, ν :: AbstractFloat = 0.0, T = Float64)


        # Define number of equations according to PDE.
        NE = PDE < 4 ? 1 : 3

        uδD = zeros(T, (N + 2, P + 1, NE))
        wδD = zeros(T, (N + 2, P + 1, NE))
        qδD = zeros(T, (N + 2, P + 1, NE))
        wBC = zeros(T, (2, P + 1, NE))
        fδD = zeros(T, (N + 2, P + 1, NE))
        rhs = zeros(T, (N + 2, P + 1, NE))
        uδI = zeros(T, (N + 1, NE))
        fδI = zeros(T, (N + 1, NE))
        qδI = zeros(T, (N + 1, NE))
        uK = zeros(T, (N + 1, 2, NE))
        fK = zeros(T, (N + 1, 2, NE))
        qK = zeros(T, (N + 1, 2, NE))
        new{T}(uδD, wδD, fδD, qδD, rhs, uδI, fδI, qδI, uK, fK, qK, wBC, T[], T[0.0], PDE, BC, γ, ν)
    end
end
"""
    init!(flow :: Flow, f :: Vector{Function}, innerPoints :: Matrix{T} where T <: AbstractFloat)

Flow initialization function. Once a `Flow` struct is created, this function initializes the `flow.uδD` values
according to a `Vector`` of analytical functions `f` and the type of equations to be solved defined in `flow.PDE`.
The number of analytical functions in `f` has to match the size of the `Flow` number of equations.
"""
function init!(flow :: Flow, f :: Vector{Function}, innerPoints :: Matrix{T} where T <: AbstractFloat)
    @assert length(f) == size(flow.wδD, 3) "Number of initial analytical functions does not match the number of primitive variables."

    for (s, fi) in enumerate(f)
        @. flow.wδD[:, :, s] = fi(innerPoints)
    end

    flow.wBC[1, :, :] .= flow.wδD[2, :, :]
    flow.wBC[2, :, :] .= flow.wδD[end - 1, :, :]

    for m ∈ 1:size(flow.uδD, 2), n ∈ 2:size(flow.uδD, 1) - 1
        u = conservative(flow.wδD[n, m, :], flow.PDE)
        flow.uδD[n, m, :] .= u
    end
end

"""
    δt!(flow :: Flow,  δx :: Vector{T} where T<: AbstractFloat,  p :: Integer, CFL :: AbstractFloat = 0.1)

Computes a suitable time step according to the fastes wave speed within the solution and the CFL condition.
"""

# @fastmath function δt!(flow :: Flow,  δx :: Vector{T} where T<: AbstractFloat,  p :: Integer, CFL :: AbstractFloat = 0.1)
#     u = 2.0 * maximum(flow.uδD)
#     dx = maximum(@. δx / (p + 1))
#     δt = CFL * dx / u # This is very restrictive, instead should be looking elementwise: Δx[i] / u[i]
#     push!(flow.δt, δt)
# end

@fastmath function δt!(flow :: Flow,  δx :: Vector{T} where T<: AbstractFloat,  p :: Integer, CFL :: AbstractFloat = 0.1)
    if flow.PDE < 4 # Linear advection or Burger (Viscous) equation
        u = maximum(flow.uδD)
        δtAdv = minimum(CFL .* δx ./ (p + 1) ./ (2.0 .* u))
        if flow.ν ≉ 0.0
            δtVis = minimum(@. 0.5 * δx * δx / flow.ν)
            push!(flow.δt, 1.0 / (1.0 / δtAdv + 1.0 / δtVis))
        else
            push!(flow.δt, δtAdv)
        end
    else # Euler equations
        pressureDensityRatio = flow.wδD[2:end - 1, :, end] ./ flow.wδD[2:end - 1, :, 1]
        pRho = maximum(pressureDensityRatio)
        soundSpeed = sqrt(flow.γ * pRho) # local element sound speed
        u = maximum(flow.wδD[2:end - 1, :, 2]) # local element velocity
        S = max(abs(u + soundSpeed), abs(u), abs(u - soundSpeed))
        δtAdv = CFL * δx[1] / (2.0 * S) # δx[1] assumes all element have same size
        push!(flow.δt, δtAdv)
    end
end

"""
    solutionBoundaryConditions!(flow :: Flow)

Applies the boundary conditions to the discontinuous conserved quantities (`flow.uδD`) according to
the `flow.BC` flag.
"""
@fastmath function solutionBoundaryConditions!(flow :: Flow)
    if flow.BC == 1 # uL 
        flow.uδD[1, :, :] .= uL_inlet
        flow.uδD[end, :, :] = flow.uδD[end - 1, :, :]
    elseif flow.BC == 2 # uG
        flow.uδD[1, :, :] .= uG_inlet
        flow.uδD[end, :, :] = flow.uδD[end - 1, :, :]
    elseif flow.BC == 3 # Liquid Holdup
        flow.uδD[1, :, :] .= αL_inlet
        flow.uδD[end, :, :] = flow.uδD[end - 1, :, :]
    elseif flow.BC == 4 # Void Fraction 
        flow.uδD[1, :, :] .= αG_inlet
        flow.uδD[end, :, :] = flow.uδD[end - 1, :, :]
    elseif flow.BC == 5 # Entalpia 
        flow.uδD[1, :, :] .= hL_intet
        flow.uδD[end, :, :] = flow.uδD[end - 1, :, :]
    else
    end
end

"""
    fluxesBoundaryConditions!(flow :: Flow)

Applies the boundary conditions to the discontinuous advection and viscous fluxes (`flow.fδD` and `flow.qδD`)
according to the `flow.BC` flag.
"""

"""
    reconstructAtFace!(flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)

Reconstruct the discontinuous conserved quantities (`flow.uδD`) at element faces, yielding `flow.uK`
"""
@fastmath function reconstructSolutionAtFace!(flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)
    @inbounds for s ∈ 1:size(flow.uδD, 3), c ∈ 2:size(flow.fδI, 1) - 1
        elementIndexL, elementIndexR = face2elements(c)
        flow.uK[c, 1, s] = @views sum(flow.uδD[elementIndexL, :, s] .* basisFace[2, :])
        flow.uK[c, 2, s] = @views sum(flow.uδD[elementIndexR, :, s] .* basisFace[1, :])
    end
    @inbounds @simd for s ∈ 1:size(flow.uδD, 3)
        flow.uK[1, 2, s] = @views sum(flow.uδD[2, :, s] .* basisFace[1, :])
        flow.uK[end, 1, s] = @views sum(flow.uδD[end - 1, :, s] .* basisFace[2, :])
    end
    flow.uK[1, 1, :] = flow.uδD[1, 1, :]
    flow.uK[1, 2, :] .= flow.uδD[1, 1, :]
    flow.uK[end, 2, :] .= flow.uK[end, 1, :]
    
end

"""
    reconstructFluxesAtFace!(flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)

Reconstruct the discontinuous fluxes (`flow.fδD` and `flow.qδD``) at element faces, yielding `flow.fK` and `flow.qK`
"""
@fastmath function reconstructFluxesAtFace!(α, ρ, α1, ρ1, u1,α2, ρ2, u2, flow :: Flow, basisFace :: Matrix{T} where T<: AbstractFloat)
    # Advection flux
    @inbounds for s ∈ 1:size(flow.uδD, 3), c ∈ 2:size(flow.fδI, 1) - 1
        elementIndexL, elementIndexR = face2elements(c)
        flow.fK[c, 1, s] = @views sum(flow.fδD[elementIndexL, :, s] .* basisFace[2, :])
        flow.fK[c, 2, s] = @views sum(flow.fδD[elementIndexR, :, s] .* basisFace[1, :])
    end
    @inbounds @simd for s ∈ 1:size(flow.uδD, 3)
        flow.fK[1, 2, s] = @views sum(flow.fδD[2, :, s] .* basisFace[1, :])
        flow.fK[end, 1, s] = @views sum(flow.fδD[end - 1, :, s] .* basisFace[2, :])
    end
    if flow.PDE == 3
        flow.fK[1, 1, :] .= primitiveEnergy(α1[1], ρ1[1], u1[1],α2[1], ρ2[1], u2[1], flow.uδD[1, 1, 1])
    else
        flow.fK[1, 1, :] .= primitive(α[1,1], ρ[1], flow.uK[1, 1, 1], flow.PDE)
    end
    flow.fK[1, 2, :] = flow.fK[1, 1, :]
    flow.fK[end, 2, :] .= flow.fK[end, 1, :]
    #
    # Diffusive flux
    if flow.ν ≉ 0.0
        @inbounds for s ∈ 1:size(flow.uδD, 3), c ∈ 2:size(flow.fδI, 1) - 1
            elementIndexL, elementIndexR = face2elements(c)
            flow.qK[c, 1, s] = @views sum(flow.qδD[elementIndexL, :, s] .* basisFace[2, :])
            flow.qK[c, 2, s] = @views sum(flow.qδD[elementIndexR, :, s] .* basisFace[1, :])
        end
        @inbounds @simd for s ∈ 1:size(flow.uδD, 3)
            flow.qK[1, 2, s] = @views sum(flow.qδD[2, :, s] .* basisFace[1, :])
            flow.qK[end, 1, s] = @views sum(flow.qδD[end - 1, :, s] .* basisFace[2, :])
        end
        if flow.BC == 0
            flow.qK[1, 1, :] .= flow.qK[end, 1, :]
            flow.qK[end, 2, :] .= flow.qK[1, 2, :]
        end
    end
end

"""
    viscousFlux!(flow :: Flow, basisFace :: Array{T, 2} where T <: AbstractFloat)

Computes the discontinuous viscous flux `flow.qδD`, where ``q = \\partial u / \\partial x``.
Note that the interaction flux `flow.uδI` is first computed, which is the solution ot the Riemann problem
at the face given the left and right solution `flow.uK`.
"""
@fastmath function viscousFlux!(flow :: Flow, dhL :: Vector{T}, dhR :: Vector{T}, Dij :: Matrix{T}, JijInv :: Vector{T}) where T <: AbstractFloat
    # Compute interaction solution
    @inbounds @simd for c ∈ 1:size(flow.uδI, 1)
        flow.uδI[c, :] .= centralFlux.(flow.uK[c, 1, :], flow.uK[c, 2, :])
    end
    # Compute qδD
    @inbounds for s ∈ 1:size(flow.uδD, 3), n ∈ 2:size(flow.uδD, 1) - 1
        faceIndexL, faceIndexR = element2faces(n)
        # Jumps of conservative quantities at faces
        jumpL = flow.uδI[faceIndexL, s] - flow.uK[faceIndexL, 2, s]
        jumpR = flow.uδI[faceIndexR, s] - flow.uK[faceIndexR, 1, s]
        # Compute RHS for reference coordinates
        flow.qδD[n, :, s] .= Dij * flow.uδD[n, :, s] + jumpL .* dhL .+ jumpR .* dhR
        # Scale dynamic viscosity, and with Jij to obtain in physical coordinates
        flow.qδD[n, :, s] .*= -JijInv[n - 1] * flow.ν
    end
end

"""
    riemannFluxes!(flow :: Flow)

Computes the interaction (aka common, numerical) advective and viscous fluxes, ie `flow.fI` and `flow.qI`,
given the left and right state of these fluxes.
Note that `fK`, where `K` stands for left `L` or right `R`, are computed as ``f(uK)``, and not using the Lagrange polynomials at x=[-1, 1].
Eg: fL = f(-1) = sum(fj * lj(-1)), as explained in Huynh original paper.
Also requires entropy fix to mitigate Gibbs oscillations (not implemented yet).
"""
@fastmath function riemannFluxes!(αØ, ρ, wavev, flow :: Flow)
    
    @inbounds @simd for c ∈ 2:size(flow.fδI, 1)-1
        uL = flow.uK[c, 1, 1]
        uR = flow.uK[c, 2, 1]
        fL = flow.fK[c, 1, 1]
        fR = flow.fK[c, 2, 1]

        uL_Post = flow.uK[c+1, 1, 1]
        uR_Ant = flow.uK[c-1, 2, 1]
        fL_Post = flow.fK[c+1, 1, 1]
        fR_Ant = flow.fK[c-1, 2, 1]


        a1 = abs((fL_Post-fR)/(uL_Post-uR+1e-10))
        a2 = abs((fR-fL)/(uR-uL+1e-10))
        a3 = abs((fL-fR_Ant)/(uL-uR_Ant+1e-10))

        
        if flow.PDE == 1 # u
            u = (uR+uL)/2
            @. flow.fδI[c, :] = montanteflux(u, fL,fR)#roeFlux(uL, uR, fL, fR, 1.0) #LLF(uL, uR, fL, fR, a1, a2, a3) #(fL+fR)/2#@. flow.fδI[c, :] = LLF(uL, uR, fL, fR, a1, a2, a3)
        
        elseif flow.PDE == 2 # 
            u = (fL+fR)/2
            @. flow.fδI[c, :] = montanteflux(u, fL,fR)# roeFlux(uL, uR, fL, fR, 1.0) #(fL+fR)/2##roeFlux(uL, uR, fL, fR, 0.0)#LLF(uL, uR, fL, fR, a1, a2, a3)
        elseif flow.PDE == 3 # 
            u = (fL+fR)/2
            @. flow.fδI[c, :] = (fL+fR)/2 - wavev*(uR-uL) #LLF
            #=elseif flow.PDE == 3 # Burgers viscous equation
            @. flow.fδI[c, :] = roeFlux(uL, uR, fL, fR, uL[1])
            @. flow.qδI[c, :] = centralFlux(qL, qR)
        #
        elseif flow.PDE == 3 # 
            @. flow.fδI[c, :] = (uL+uR)/2#LLF(uL, uR, fL, fR, a1, a2, a3)
        elseif flow.PDE == 4 # 
            @. flow.fδI[c, :] = (uL+uR)/2#LLF(uL, uR, fL, fR, a1, a2, a3)
            
        else # Euler equations
            # flow.fδI[c, :] .= rusanov(uL, uR, fL, fR, wL, wR, flow.γ)
            # flow.fδI[c, :] .= hll(uL, uR, fL, fR, wL, wR, flow.γ)
            flow.fδI[c, :] .= hllc(uL, uR, fL, fR, wL, wR, flow.γ)
        =#
        end
    end
    flow.fδI[1] = flow.fδD[1, 1,1]
    flow.fδI[end] = flow.fδI[end-1]
end

"""
    rhs!(flow :: Flow, mesh :: Mesh)

Compute the RHS of the equations so that the solution can be advanced using a time integration methods such as Euler or Runge-Kutta.
The multiple steps are defined next:
1. Apply the boundary conditions to the discontinuous solution `flow.uδD`.
2. Reconstruct discontinuous solution `flow.uδD` at faces, obtaining `flow.uK`.
3. Compute the viscous flux `flow.qδD`.
4. Compute the advection flux `flow.fδD`.
5. Reconstruct advection and viscous fluxes at faces, obtaining `flow.fK` and `flow.qK`.
6. Compute common fluxes at element faces `flow.fδI` and `flow.qδI`.
7. Compute the equation RHS, du/dt = RHS = -1/Jij * df/dr, where `Jij` is the Jacobian of the mesh.
"""
@fastmath function rhs!(flow :: Flow, αØ, ρ, mesh :: Mesh, α0, ρ0, Δx, p, θ, dif_art, veloc,fi0,τi,τw,αGF0,D,a1,a2,a3,a4,a5,a6,a7,a8,α1, ρ1, u1,α2, ρ2, u2,wavev, flowpat, Sw, Si)
    # Get correction function derivative and the Lagrange derivative matrix at element integration points.
    basisFace = mesh.standardElement.basisFace
    dhL = mesh.polynomialData.correctionDerivativeL
    dhR = mesh.polynomialData.correctionDerivativeR
    Dij = mesh.polynomialData.Dij
    JijInv = mesh.JijInv

    # 1. Apply BCs to solution: flow.uδD
    solutionBoundaryConditions!(flow)

    # 2. Extrapolate conserved quantities at faces: flow.uK
    reconstructSolutionAtFace!(flow, basisFace)


    # 3. Compute advection flux at solution points: flow.fδD
    velocidade = 0
    @inbounds for n ∈ 2:size(flow.uδD, 1)-1
        for m ∈ 1:size(flow.uδD, 2) 
            #=
            posicao = Δx*(n-2) + (m-1)*Δx
            # Velocidades Analíticas
            if posicao <= (10*t+0.5*9.81*t^2)
                if flow.BC == 3
                    velocidade = (2*9.81*posicao+100)^(1/2)
                elseif flow.BC == 4
                    velocidade = 0
                end
            else
                if flow.BC == 3
                    velocidade = 10 + 9.81*t
                elseif flow.BC == 4
                    velocidade = -(1-0.2)*9.81*t/0.2
                end
            end
            =#
            # Argumento de velocidade só funciona pra Lobatto Points
            if flow.PDE == 3
                flow.fδD[n, m, 1] = primitiveEnergy(α1[n-1], ρ1[n-1], u1[n-1],α2[n-1], ρ2[n-1], u2[n-1], flow.uδD[n, m, 1])
            else
                flow.fδD[n, m, 1] = primitive(veloc[n,m], ρ[n-1], flow.uδD[n, m, 1], flow.PDE) #primitive(α[n, m, 1], ρ[aux], flow.uδD[n, m, 1], flow.PDE)
            end
        end
    end
    flow.fδD[1, :, 1] = flow.fδD[2, :, 1]
    flow.fδD[size(flow.uδD, 1), :, 1] = flow.fδD[size(flow.uδD, 1)-1, :, 1]

    

    # 4. Extrapolate advection and diffusion fluxes at faces: flow.fK, flow.qK
    reconstructFluxesAtFace!(αØ, ρ, α1, ρ1, u1,α2, ρ2, u2, flow, basisFace)

    

    # 5. Compute common advection and diffusion fluxes at faces: flow.fδI, flow.qδI
    riemannFluxes!(αØ, ρ, wavev, flow)

    

    # 6. Compute the equation RHS
    @inbounds for s ∈ 1:size(flow.uδD, 3), n ∈ 2:size(flow.uδD, 1) - 1
        faceIndexL, faceIndexR = element2faces(n)
        # Advection flux
        jumpLAdv = flow.fδI[faceIndexL, s] - flow.fK[faceIndexL, 2, s]
        jumpRAdv = flow.fδI[faceIndexR, s] - flow.fK[faceIndexR, 1, s]
        # Compute RHS for reference coordinates
        flow.rhs[n, :, s] .= Dij * flow.fδD[n, :, s] + jumpLAdv .* dhL .+ jumpRAdv .* dhR
        # Scale with Jij to obtain in physical coordinates
        flow.rhs[n, :, s] .*= -JijInv[n - 1]
        
        
    end
    # Adicionando demais termos (eu que fiz):
        #
        tamanho = size(flow.uδD[:,1,1],1)
        #
        if flow.PDE == 1
            # Gradientes de pressão (eu que fiz) = 
            global gradp = zeros(tamanho-2,1)
            for i=2:(tamanho-3)
                gradp[i] = (p[i]-p[i-1])/Δx
            end
            gradp[tamanho-2] = gradp[tamanho-3]
            #gradp[tamanho-1] = gradp[tamanho-3]

            # Adicionando demais termos (eu que fiz):
            if flowpat == 1
                if flow.BC == 1
                    flow.rhs[2:(end-1), :, 1] = ρ0[:].*α0[2:end - 1,:].*flow.uδD[2:(end-1), :, 1].*flow.rhs[2:(end-1), :, 1] .- α0[2:end - 1,:].*gradp[:] .+ α0[2:end - 1,:].*ρ0[:]*(-9.81)*sin(θ) .- dif_art[:] .+ 2*τi[:].*(αGF0[:]).^(1/2)./(D/2) .- (4/D)*τw[:]
                else
                    flow.rhs[2:(end-1), :, 1] = ρ0[:].*α0[2:end - 1,:].*flow.uδD[2:(end-1), :, 1].*flow.rhs[2:(end-1), :, 1] .- α0[2:end - 1,:].*gradp[:] .+ α0[2:end - 1,:].*ρ0[:]*(-9.81)*sin(θ) .- dif_art[:] .- 2*αGF0[:].*τi[:]./((αGF0[:]).^(1/2).*D/2)
                end
            elseif flowpat == 2
                if flow.BC == 1
                    flow.rhs[2:(end-1), :, 1] = ρ0[:].*α0[2:end - 1,:].*flow.uδD[2:(end-1), :, 1].*flow.rhs[2:(end-1), :, 1] .- α0[2:end - 1,:].*gradp[:] .+ α0[2:end - 1,:].*ρ0[:]*(-9.81)*sin(θ) .- dif_art[:] .+ 2*τi[:].*Si[:]./(3.14*D^2/4) .- τw[:].*Sw[:]./(3.14*D^2/4)
                else
                    flow.rhs[2:(end-1), :, 1] = ρ0[:].*α0[2:end - 1,:].*flow.uδD[2:(end-1), :, 1].*flow.rhs[2:(end-1), :, 1] .- α0[2:end - 1,:].*gradp[:] .+ α0[2:end - 1,:].*ρ0[:]*(-9.81)*sin(θ) .- dif_art[:] .- 2*τi[:].*Si[:]./(3.14*D^2/4) .- τw[:].*Sw[:]./(3.14*D^2/4)
                end 
            end  
        elseif flow.PDE == 3
            flow.rhs[2:(end-1), :, 1] = flow.rhs[2:(end-1), :, 1] .+ a4[:] .+ a5[:] .+ a6[:] .+ a7[:]

        end
    
    # 7. Extrapolate conserved quantities at faces: flow.uK
    reconstructSolutionAtFace!(flow, basisFace)

    # 8. Pegando os argumentos nas faces (eu que fiz)
    for c ∈ 1:size(flow.uδI, 1)
        flow.uδI[c, :] .= centralFlux.(flow.uK[c, 1, :], flow.uK[c, 2, :])
    end

end

"""
    stepRK3SSPS3!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)

Update the numerical solution `flow.uδD` using the Runge-Kutta strong stability preserving (SSP) 3rd order (3 stages) scheme.
The time step is computed using the CFL condition `CFL :: AbstractFloat`.
"""
@fastmath function stepRK3SSPS3!(flow :: Flow, u0, α, ρ, α0, ρ0, Δt,Δx,P,θ,dif_art,vel,Ø, mesh :: Mesh, CFL :: AbstractFloat = 0.1)
    N = length(ρ0)
    #Δt = Δt/2

    # Termo N (paper de Eduarda) inicial:
    if flow.PDE == 1 
      
        u0RK = α0[:].*ρ0[:].*u0[2:end - 1,:] 

        divisor = α0[:].*ρ0[:]
   elseif flow.PDE == 2
        u0RK = u0[2:end - 1,:] 
        divisor = 1 
   end

    
    rhs!(flow, α, ρ, mesh, α0, ρ0, Δx, P, θ,dif_art,vel)



    u1RK = @. u0RK + Δt*flow.rhs[2:end - 1, :, :]
   
    u1RK[end,:] = u1RK[end-1,:]
    

    
    @. flow.uδD[2:end - 1, :, :] = u1RK ./ divisor

    #Limitação:
    uc = flow.uδD[2:end - 1, :, 1]'
    flow.uδD[2:end - 1, :, 1] = hMLP(uc, Grau+1, N, V, Grau)


    # Condição de Contorno Novamente:
    flow.uδD[2, :, :] = flow.uδD[1, :, :]
    flow.uδD[end, :, :] = flow.uδD[end-1, :, :]


    rhs!(flow, α, ρ, mesh, α0, ρ0,Δx,P,θ,dif_art,vel)

    u2RK = @. 0.75 * u0RK + 0.25 * (u1RK + Δt * flow.rhs[2:end - 1, :, :])

    u2RK[end,:] = u2RK[end-1,:]

    @. flow.uδD[2:end - 1, :, :] = u2RK ./ divisor

    #Limitação:
    uc = flow.uδD[2:end - 1, :, 1]'
    flow.uδD[2:end - 1, :, 1] = hMLP(uc, Grau+1, N, V, Grau)


    # Condição de Contorno Novamente:
    flow.uδD[2, :, :] = flow.uδD[1, :, :]
    flow.uδD[end, :, :] = flow.uδD[end-1, :, :]

 

    rhs!(flow, α, ρ, mesh, α0, ρ0,Δx,P,θ,dif_art,vel)

    @. flow.uδD[2:end - 1, :, :] = 1.0 / 3.0 * u0RK + 2.0 / 3.0 * (u2RK  + Δt * flow.rhs[2:end - 1, :, :])
    flow.uδD[end - 1, :, :] = flow.uδD[end - 2, :, :]
    @. flow.uδD[2:end - 1, :, :] = flow.uδD[2:end - 1, :, :] ./ divisor

    #Limitação:
    uc = flow.uδD[2:end - 1, :, 1]'
    flow.uδD[2:end - 1, :, 1] = hMLP(uc, Grau+1, N, V, Grau)

    # Condição de Contorno Novamente:
    flow.uδD[2, :, :] = flow.uδD[1, :, :]
    flow.uδD[end, :, :] = flow.uδD[end-1, :, :]



end

"""
    stepEuler!(flow :: Flow, mesh :: Mesh, CFL :: AbstractFloat = 0.1)

Update the numerical solution `flow.uδD` using the Euler scheme.
The time step is computed using the CFL condition `CFL :: AbstractFloat`.
"""
@fastmath function stepEuler!(flow :: Flow, u0, α, ρ, α0, ρ0, Δt,Δx,P,θ,dif_art,vel,Ø,fi0,τi,τw,αGF0,D,a1,a2,a3,a4,a5,a6,a7,a8,α1, ρ1, u1,α2, ρ2, u2,wavev, flowpat, Sw, Si, mesh :: Mesh, CFL :: AbstractFloat = 0.1)
    N = length(ρ0)
    #Δt = Δt/2

    # Termo N (paper de Eduarda) inicial:
    if flow.PDE == 1 
        u0RK = α0[2:end - 1,:].*ρ0[:].*u0[2:end - 1,:] 
        divisor = α0[2:end - 1,:].*ρ0[:]

   elseif flow.PDE == 2
        u0RK = u0[2:end - 1,:] 
        divisor = 1 

    elseif flow.PDE == 3
        u0RK = a1[:] + a2[:] 
        divisor = a8[:] 
   end

   
    rhs!(flow, α, ρ, mesh, α0, ρ0, Δx, P, θ,dif_art,vel,fi0,τi,τw,αGF0,D, a1,a2,a3,a4,a5,a6,a7,a8, α1, ρ1, u1,α2, ρ2, u2,wavev, flowpat, Sw, Si)

    

    u1RK = @. u0RK .+ Δt*flow.rhs[2:end - 1, :, :]

   
    u1RK[end,:] = u1RK[end-1,:]

    
    @. flow.uδD[2:end - 1, :, :] = u1RK ./ divisor

    #Limitação:
    uc = flow.uδD[2:end - 1, :, 1]'
    if flow.PDE == 3
        flow.uδD[2:end - 1, :, 1] = hMLP(uc, Grau+1, N, V, Grau)
    else
        flow.uδD[2:end - 1, :, 1] = hmlp2(uc, Grau+1, N, V, Grau) # 
    end
    
    # Condição de Contorno Novamente:
    flow.uδD[2, :, :] = flow.uδD[1, :, :]
    flow.uδD[end, :, :] = flow.uδD[end-1, :, :]

    
end

"""
    flattenSolution(u :: Array{T, 2} where T <: AbstractFloat)

Returns a `Vector` type containing flow.uδD with all integration points of the mesh flattened in a single `Vector`.
"""
@fastmath function flattenSolution(u :: Array{T, 2} where T <: AbstractFloat)
    uFlat = zeros(eltype(u), (size(u, 1) - 2) * size(u, 2))
    k = size(u, 2) # number of solution points per element
    @inbounds  for n ∈ 2:size(u, 1) - 1
        uFlat[(n - 2) * k + 1:(n - 2) * k + k] = u[n, :]
    end
    return uFlat
end
#
"""
    function updatePrimitiveVars!(flow :: Flow)

Update the primitive variables array of the a `Flow`, ie `Flow.wδD`.
"""
@fastmath function updatePrimitiveVars!(flow :: Flow)
    @inbounds for m ∈ 1:size(flow.uδD, 2), n ∈ 2:size(flow.uδD, 1) - 1
        flow.wδD[n, m, :] .= primitive(flow.uδD[n, m, :], flow.PDE)
    end
end
