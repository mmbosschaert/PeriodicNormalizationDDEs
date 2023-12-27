function locate_double_hopf(branch)

    num_unstable = zeros(Int,length(branch))
    for (i,p) in enumerate(branch)
        p.stability
        # remove ω from stability field
        ind1 = findfirst(λ -> imag(λ) ≈ p.ω, p.stability)
        ind2 = findfirst(λ -> imag(λ) ≈ -p.ω, p.stability)
        eigenvalues = p.stability[1:end .!= ind1 .&& 1:end .!= ind2]
        # count unstable eigenvalues
        num_unstable[i] = sum(real(eigenvalues) .> 0.0)
    end
    ind_double_hopf = findall(p -> abs(p) == 2, num_unstable[2:end] - num_unstable[1:end-1])
    ind_double_hopf
end

function normalform(jet, hoho::DoubleHopf, τs)

    m = length(τs)
    φ = repeat(hoho.coords,1,m)
    α = hoho.parameters

    λ₁ = hoho.ω₁*im
    λ₂ = hoho.ω₂*im

    Δ(λ) = jet.Δ(φ,α,λ)
    Δ′(λ)= jet.Δ′(φ,α,λ)
    # q = nullspace(Δ(λ); atol=1e-14)
    # p = nullspace(Δ(λ)';atol=1e-14)

    _, s, V = svd(Δ(λ₁))
    indxmin = last(findmin(s))
    q1 = V[:, indxmin]

    _, s, V = svd(transpose(Δ(λ₁)))
    indxmin = last(findmin(s))
    p1 = V[:, indxmin]

    _, s, V = svd(Δ(λ₂))
    indxmin = last(findmin(s))
    q2 = V[:, indxmin]

    _, s, V = svd(transpose(Δ(λ₂)))
    indxmin = last(findmin(s))
    p2 = V[:, indxmin]

    # normalize
     
    p1 /=  transpose(p1)*(Δ′(λ₁)*q1)
    p2 /=  transpose(p2)*(Δ′(λ₂)*q2)

    q = [q1, q2]
    p = [p1, p2]
    # first(transpose(p1)*(Δ′(λ₁)*q1)) ≈ 1.0
    # first(transpose(p2)*(Δ′(λ₂)*q2)) ≈ 1.0

    # multi-linear forms at bt point
    Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
    B(v₁,v₂)  = jet.D2(φ,α,Ξ(v₁),Ξ(v₂))
    C(v₁,v₂,v₃)  = jet.D3(φ,α,Ξ(v₁),Ξ(v₂),Ξ(v₃))
    A1(v₁,p₁)  = jet.D11(φ,α,Ξ(v₁),p₁)
    J1  = jet.D01(φ,α)

    ϕ1(θ) = exp(λ₁*θ)*q1
    ϕ2(θ) = exp(λ₂*θ)*q2
    ϕs = [ϕ1, ϕ2]

    # Quadratic center manifold
    h1100(_) = Δ(0)\B(ϕ1,conj∘ϕ1)
    h2000(θ) = exp(2*λ₁*θ)*(Δ(2*λ₁)\B(ϕ1,ϕ1))
    h1010(θ) = exp((λ₁+λ₂)*θ)*(Δ(λ₁+λ₂)\B(ϕ1,ϕ2))
    h1001(θ) = exp((λ₁-λ₂)*θ)*(Δ(λ₁-λ₂)\B(ϕ1,conj∘ϕ2))
    h0020(θ) = exp(2λ₂*θ)*(Δ(2λ₂)\B(ϕ2,ϕ2))
    h0011(_) = Δ(0)\B(ϕ2,conj∘ϕ2)

    g2100 = (1/2)*transpose(p1)*( 2*B(h1100,ϕ1) + B(h2000, conj∘ϕ1) + C(ϕ1,ϕ1,conj∘ϕ1) )
    g1011 = transpose(p1)*( B(h0011,ϕ1) + B(h1001,ϕ2) + B(h1010,conj∘ϕ2) + C(ϕ1,ϕ2,conj∘ϕ2) )
    g1110 = transpose(p2)*( B(conj∘h1001,ϕ1) + B(h1010,conj∘ϕ1)+ B(h1100,ϕ2) + C(ϕ1,conj∘ϕ1,ϕ2) )
    g0021 = (1/2)*transpose(p2)*( 2*B(h0011,ϕ2) + B(h0020,conj∘ϕ2) + C(ϕ2,ϕ2,conj∘ϕ2) )

    θ = real(g1011)/real(g0021)
    δ = real(g1110)/real(g2100)

    Id = diagm(ones(2))
    Γ = [transpose(p[i])*(A1(ϕs[i],Id[:,j])+B(ϕs[i],_ -> Δ(0)\(J1*Id[:,j]))) for i=1:2, j=1:2]
    K = inv(real(Γ))

    h000001 = Δ(0)\(J1*K[:,1])
    h000010 = Δ(0)\(J1*K[:,2])
    h0000 = [_ -> h000001, _ -> h000010]
    b = [imag(transpose(p[i]) * (A1(ϕs[i],K[:,j]) + B(ϕs[i],h0000[j]))) for i=1:2, j=1:2]

    nmfm = DoubleHopfNormalform(
        q,
        g2100,
        g1011,
        g1110,
        g0021,
        θ,
        δ,
        b,
        K,
        h0000,
        h0011,
        h0020,
        h2000,
        h1100
        )

    hoho = @set hoho.nmfm = nmfm
end


