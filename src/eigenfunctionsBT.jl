function eigenfunctionsBT(Δ, Δs)

    Δ′, Δ′′, Δ′′′, Δ⁴, Δ⁵ = Δs 

    # construct null Jordan chains
    q₀ = vec(nullspace(Δ(0)))
    p₁ = vec(nullspace(Δ(0)'))
    q₁ = -([Δ(0) p₁; [q₀' 0]]\[Δ′*q₀; 0])[1:end-1]
    p₀ = -([Δ(0) p₁; [q₀' 0]]'\[Transpose(Δ′)*p₁; 0])[1:end-1]
    p₁ = p₁'
    p₀ = p₀'

    n = length(q₀)

    @assert ≈(Δ(0)*q₀, zeros(n,1); atol=1e-14)
    @assert ≈(Δ′*q₀ + Δ(0)*q₁, zeros(n,1); atol=1e-14)
    @assert ≈(p₁*Δ(0), zeros(1,n) ; atol=1e-14)
    @assert ≈(p₁*Δ′ + p₀*Δ(0), zeros(1,n); atol=1e-14)

    # normalize eigenfunctions
    c₁ = p₀*Δ′*q₀ + 1/2*p₁*Δ′′*q₀
    p₀ /= c₁
    p₁ /= c₁
    c₂ = p₀*Δ′*q₁ + 1/2*p₀*Δ′′*q₀ + 1/2*p₁*Δ′′*q₁ + 1/6*p₁*Δ′′′*q₀
    p₀ -= c₂*p₁
    @assert ≈(p₀*Δ′*q₀ + 1/2*p₁*Δ′′*q₀, 1)
    @assert ≈((p₀*Δ′*q₁ + 1/2*p₀*Δ′′*q₀ + 1/2*p₁*Δ′′*q₁ + 1/6*p₁*Δ′′′*q₀), 0.0; atol=1e-14)

    # define eigenfucntions
    φ₀ = Polynomial([q₀], :ϑ)
    φ₁ = Polynomial([q₁,q₀], :ϑ)

    # define pairings with left eigenfucntions
    ψ₁(φ::Polynomial) = p₁*sum(map((i,δ,c) -> δ*c/i, 1:length(φ), Δs, φ.coeffs))
    ψ₀(φ::Polynomial) = p₀*sum(map((i,δ,c) -> δ*c/i, 1:length(φ), Δs, φ.coeffs)) + 
                            p₁*sum(map((i,δ,c) -> δ*c/(i*(i+1)), 1:length(φ), Δs[2:end], φ.coeffs))

    # assert normalized eigenfucntions pairing
    @assert ≈(ψ₁(φ₁), 1.0; atol=1e-14)
    @assert ≈(ψ₀(φ₀), 1.0; atol=1e-14)
    @assert ≈(ψ₁(φ₀), 0.0; atol=1e-14)
    @assert ≈(ψ₀(φ₁), 0.0; atol=1e-14)

    p₁, p₀, q₀, q₁, ψ₁, ψ₀, φ₀, φ₁

end
