mutable struct Hopf
    coords::Vector{Float64}
    parameters::Vector{Float64}
    v::Union{Vector{ComplexF64},Vector{Complex{Num}}}
    ω::Float64
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{ComplexF64,Nothing}
end
Hopf(coords,parameters,v,ω) = Hopf(coords,parameters,v,ω,nothing,nothing)

mutable struct DoubleHopfNormalform
    q
    g2100
    g1011
    g1110
    g0021
    θ
    δ
    b
    K
    h0000
    h0011
    h0020
    h2000
end

mutable struct DoubleHopf
    coords::Vector{Float64}
    parameters::Vector{Float64}
    v::Union{Vector{ComplexF64},Vector{Complex{Num}}}
    ω₁::Float64
    ω₂::Float64
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{DoubleHopfNormalform,Nothing}
end
DoubleHopf(coords,parameters,v,ω₀,ω₁) = DoubleHopf(coords,parameters,v,ω₀,ω₁,nothing,nothing)

function piont_to_hoho(p)
    if p.stability === nothing
        println("Need to calculate stability")
        return nothing
    else
        ω_indx = findall( x -> abs(real(x)) < 1e-03, p.stability)
        if length(ω_indx) !== 4
            println("Unable to extract frequencies")
        end
        ω₀, ω₁ = sort(abs.(imag(p.stability[ω_indx])))[[1,3]]
        DoubleHopf(p.coords,p.parameters,[0.0+0.0im],ω₀,ω₁)
    end
end

function vec(hopf::Hopf,_)::Vector{Float64}
    vcat(hopf.coords, hopf.parameters,real(hopf.v),imag(hopf.v),hopf.ω)
end

function vec_to_point(v::Vector{Float64},::Hopf,_)
    dims = div(length(v) - 1,4)
    Hopf(v[1:dims], v[dims+1:2dims], v[2dims+1:3dims] + v[3dims+1:4dims]*im, v[4dims+1])
end

function Hopf_res!(res, model, τs, Δre, Δim, xx, xx_prev, n)

    x, α, vre, vim , ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2], xx[2n+3:3n+2], xx[3n+3]
    vre_prev, vim_prev = xx_prev[n+3:2n+2], xx_prev[2n+3:3n+2]

    res[1:n] = model(repeat(x,1,length(τs)),α)
    res[n+1:2n]  = Δre(repeat(x,1,length(τs)),α,ω)*vre - Δim(repeat(x,1,length(τs)),α,ω)*vim
    res[2n+1:3n] = Δre(repeat(x,1,length(τs)),α,ω)*vim + Δim(repeat(x,1,length(τs)),α,ω)*vre
    res[3n+1] = dot(vre,vre_prev) - dot(vim,vim_prev) - 1
    res[3n+2] = dot(vim,vre_prev) + dot(vre,vim_prev)
end

function Hopf_res!(res, model, τs, Δre, Δim, xx, hopf_prev::Hopf, n)
    x, α, vre, vim , ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2], xx[2n+3:3n+2], xx[3n+3]
    vre_prev, vim_prev = real(hopf_prev.v), imag(hopf_prev.v)

    res[1:n] = model(repeat(x,1,length(τs)),α)
    res[n+1:2n]  = Δre(repeat(x,1,length(τs)),α,ω)*vre - Δim(repeat(x,1,length(τs)),α,ω)*vim
    res[2n+1:3n] = Δre(repeat(x,1,length(τs)),α,ω)*vim + Δim(repeat(x,1,length(τs)),α,ω)*vre
    # res[3n+1] = dot(vre,vre_prev) - dot(vim,vim_prev) - 1
    res[3n+1] = dot(vre,vre) + dot(vim,vim) - 1
    res[3n+2] = dot(vim,vre_prev) + dot(vre,vim_prev)
end


function Hopf_res(model, τs, Δre, Δim, xx, hopf_prev::Hopf, n)
    res = zeros(3n+3)
    Hopf_res!(res, model, τs, Δre, Δim, xx, hopf_prev::Hopf, n)
    res
end

# function Hopf_res!(res, model, τs, Δ, xx, xx_prev, n)
# 
#     x, α, v, ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2] + xx[2n+3:3n+2]*im, xx[3n+3]
#     v_prev = xx_prev[n+3:2n+2] + xx_prev[2n+3:3n+2]*im
# 
#     res[1:n] = model(repeat(x,1,length(τs)),α)
#     res[n+1:2n]  = real(Δ(repeat(x,1,length(τs)),α,ω*im)*v)
#     res[2n+1:3n] = imag(Δ(repeat(x,1,length(τs)),α,ω*im)*v)
#     res[3n+1] = real(dot(v,v_prev)) - 1
#     res[3n+2] = imag(dot(v,v_prev))
# end

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
    # @show p1
    # @show p2
    # @show transpose(p1)*(Δ(λ₁))
    # @show transpose(p2)*(Δ(λ₂))
    # @show first(transpose(p1)*(Δ′(λ₁)*q1)) ≈ 1.0
    # @show first(transpose(p2)*(Δ′(λ₂)*q2)) ≈ 1.0

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
    # h0000 = [Δ(0)\(J1*K[:,i]) for i=1:2]
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
        h2000)

    hoho = @set hoho.nmfm = nmfm
end

function normalform(jet, hopf::Hopf, τs)

    m = length(τs)
    φ = repeat(hopf.coords,1,m)
    α = hopf.parameters

    λ = hopf.ω*im
    Δ(λ) = jet.Δ(φ,α,λ)
    Δ′(λ)= jet.Δ′(φ,α,λ)

    q = qr(Δ(λ)').Q[:,end]
    p = qr(Δ(λ)).Q[:,end]

    _, s, V = svd(Δ(λ))
    indxmin = last(findmin(s))
    q = V[:, indxmin]

    _, s, V = svd(transpose(Δ(λ)))
    indxmin = last(findmin(s))
    p = V[:, indxmin]

    # normalize
    p /=  transpose(p)*(Δ′(λ)*q)

    # multi-linear forms at bt point
    Ξ(h) = vcat([h(-τ) for τ ∈ τs]...)
    B(v₁,v₂)  = jet.D2(φ,α,Ξ(v₁),Ξ(v₂))
    C(v₁,v₂,v₃)  = jet.D3(φ,α,Ξ(v₁),Ξ(v₂),Ξ(v₃))

    ϕ(θ) = exp(λ*θ)*q
    h2000(θ) = exp(2*λ*θ)*(Δ(2*λ)\B(ϕ,ϕ))
    h1100(θ) = Δ(0)\B(ϕ,conj∘ϕ)

    c₁ = first((1/2)*transpose(p)*(B(conj∘ϕ,h2000) + 2*B(ϕ,h1100) + C(ϕ,ϕ,conj∘ϕ)))

    ℓ₁ = real(c₁)/hopf.ω
    ℓ₁
end

function locate_genh(jet, branch, τs)
    for hopf in branch
        hopf.nmfm = normalform(jet, hopf, τs)
    end
    nmfm = [hopf.nmfm for hopf in branch]
    findall(!iszero,sign.(nmfm)[2:end] - sign.(nmfm)[1:end-1])
end

function SetupHopfBranch(jet,hopf_point,τs; parameterbounds=nothing, δ=.001, δmin=1e-06, δmax=0.01, MaxNumberofSteps = 250)
    dims = length(hopf_point.coords)

    # create f and jac for continuation of Hopf points
    @variables xx[1:3dims+3] v_prev[1:2dims]
    xx = Symbolics.scalarize(xx)
    v_prev = Symbolics.scalarize(v_prev)
    res = similar(xx)
    hopf_prev = Hopf(zeros(dims), zeros(2), v_prev[1:dims] + v_prev[dims+1:2dims]*im, 0.0, nothing, nothing)
    Δre, Δim = characteristic_matrices_unevaluated_re_im(jet.system, dims, τs)
    Hopf_res!(res, jet.system, τs, Δre, Δim, xx, hopf_prev, dims)
    res[end] = 0.0

    hopfJac = Symbolics.jacobian(res,xx)

    f  = build_function(res, xx, v_prev, expression=Val{false})[1]
    df = build_function(hopfJac,xx,v_prev, expression=Val{false})[1]

    # hopf branch I
    x₀ = vec(hopf_point,nothing)
    x₀prev = x₀
    jac = df(x₀,x₀prev)
    # tangent vector
    v₀ = qr(jac').Q[:,end]

    # @show display(f(x₀,[real(hopf_point.v); imag(hopf_point.v)]))


    # hopf_branch = (points = Hopf[], tangents = [], stepsizes = [], 
    #     f = (x,hopf) ->  f(x,[real(hopf.v);imag(hopf.v)]),
    #    df = (x,hopf) -> df(x,[real(hopf.v);imag(hopf.v)]),
    #     parameterbounds=parameterbounds,
    #     δ=δ,
    #     δmin=δmin,
    #     δmax=δmax,
    #     MaxNumberofSteps = MaxNumberofSteps,
    #     con_par = nothing 
    # )

    hopf_branch = (points = Hopf[], tangents = [], stepsizes = [], 
        f = (x,hopf) ->  Hopf_res(jet.system, τs, Δre, Δim, x, hopf, dims),
       df = (x,hopf) -> df(x,[real(hopf.v);imag(hopf.v)]),
        parameterbounds=parameterbounds,
        δ=δ,
        δmin=δmin,
        δmax=δmax,
        MaxNumberofSteps = MaxNumberofSteps,
        con_par = nothing 
    )

    push!(hopf_branch.points, hopf_point)
    push!(hopf_branch.tangents, v₀)
    push!(hopf_branch.stepsizes, 0.0)
    hopf_branch
end

function hopf_from_hoho(jet, hoho1, τs; ϵ = 1e-03)
    # TODO: add second curve
    β₁ = -ϵ
    # coords = hoho1.coords + [sqrt(-β₁/real(hoho1.nmfm.g2100)); 0.0]
    # β₂ = 0.00001
    # coords = [0.0; sqrt(-β₂/real(hoho1.nmfm.g0021))]
    coords = [0.0; 0.0]

    # hopf1 = Hopf( coords, hoho1.parameters + hoho1.nmfm.K*[0.0; β₂], hoho1.nmfm.q[2], hoho1.ω₂)
    hopf1 = Hopf(coords, hoho1.parameters + hoho1.nmfm.K*[β₁; 0.0], hoho1.nmfm.q[1], hoho1.ω₁)

    m = length(τs)
    φ = repeat(hopf1.coords,1,m)
    α = hopf1.parameters
    λ₁ = hoho1.ω₁*im
    Δ(λ) = jet.Δ(φ,α,λ)
    _, s, V = svd(Δ(λ₁))
    indxmin = last(findmin(s))
    q1 = V[:, indxmin]

    hopf_prev = hopf1
    c = sqrt(1/(dot(real(hopf1.v),real(hopf_prev.v)) - dot(imag(hopf1.v),imag(hopf_prev.v))))
    c = sqrt(1/(dot(real(hopf1.v),real(hopf1.v)) - dot(imag(hopf1.v),imag(hopf1.v))))
    hopf1 = @set hopf1.v = c*q1
end
