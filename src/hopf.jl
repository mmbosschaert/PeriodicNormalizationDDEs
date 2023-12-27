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
    h1100
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

mutable struct GenHopfNormalform
    b12
    c1
    c2
    phi
    K01
    h10mu
    h00mu
    h1100
    h2000
    h3000
    h2100
end

mutable struct GenHopf
    coords::Vector{Float64}
    parameters::Vector{Float64}
    v::Union{Vector{ComplexF64},Vector{Complex{Num}}}
    ω::Float64
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{GenHopfNormalform,Nothing}
end
GenHopf(coords,parameters,v,ω) = GenHopf(coords,parameters,v,ω,nothing,nothing)

function point_to_genhopf(p::Hopf)
    if p.stability === nothing
        println("Need to calculate stability")
        return nothing
    else
        GenHopf(p.coords,p.parameters,p.stability,p.ω)
    end
end

function point_to_hoho(p)
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

# Function to correct Hopf point in one parameter 
function Hopf_res!(res, model, τs, Δre, Δim, xx, hopf_prev::Hopf, n)
    x, α, vre, vim , ω = xx[1:n], xx[n+1:n+2], xx[n+3:2n+2], xx[2n+3:3n+2], xx[3n+3]
    vre_prev, vim_prev = real(hopf_prev.v), imag(hopf_prev.v)


    res[1:n] = model(repeat(x,1,length(τs)),α)
    res[n+1:2n]  = Δre(repeat(x,1,length(τs)),α,ω)*vre - Δim(repeat(x,1,length(τs)),α,ω)*vim
    res[2n+1:3n] = Δre(repeat(x,1,length(τs)),α,ω)*vim + Δim(repeat(x,1,length(τs)),α,ω)*vre
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
    h1100(_) = Δ(0)\B(ϕ,conj∘ϕ)

    c₁ = first((1/2)*transpose(p)*(B(conj∘ϕ,h2000) + 2*B(ϕ,h1100) + C(ϕ,ϕ,conj∘ϕ)))

    ℓ₁ = real(c₁)/hopf.ω
    ℓ₁
end

function detect_genh(jet, branch, τs)
    for hopf in branch
        hopf.nmfm = normalform(jet, hopf, τs)
    end
    nmfm = [hopf.nmfm for hopf in branch]
    findall(!iszero,sign.(nmfm)[2:end] - sign.(nmfm)[1:end-1])
end

# add to Hopf points
+(hopf1::Hopf, hopf2::Hopf) = Hopf(hopf1.coords + hopf2.coords, hopf1.parameters + hopf2.parameters, hopf1.v + hopf2.v, hopf1.ω + hopf2.ω)

# multiply Hopf point with scalar
*(a::Float64, hopf::Hopf) = Hopf(a*hopf.coords, a*hopf.parameters, a*hopf.v, a*hopf.ω)

# function to locate generalized Hopf points
# input are two Hopf points
# output is the Hopf point in between
function locate_genh(jet, hopf1, hopf2, τs)
    # calulate point in between
    hopf = 0.5*(hopf1 + hopf2)
    # correct point
    hopf = correct_hopf(jet, hopf, τs)
    # calculate normal form
    hopf.nmfm = normalform(jet, hopf, τs)
end

function defining_system_Hopf(jet, τs, dims)
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

    f, df, Δre, Δim
end

function SetupHopfBranch(jet,hopf_point,τs; parameterbounds=nothing, δ=.001, δmin=1e-06, δmax=0.01, MaxNumberofSteps = 250, NumberOfFails = 4)
    dims = length(hopf_point.coords)

    # create f, df, Δre, Δim for continuation of Hopf points
    f, df, Δre, Δim = defining_system_Hopf(jet, τs, dims)

    # hopf branch I
    x₀ = vec(hopf_point,nothing)
    x₀prev = x₀
    jac = df(x₀,x₀prev)
    # tangent vector
    v₀ = qr(jac').Q[:,end]

    hopf_branch = (points = Hopf[], tangents = [], stepsizes = [], 
        f = (x,hopf) ->  Hopf_res(jet.system, τs, Δre, Δim, x, hopf, dims),
       df = (x,hopf) -> df(x,[real(hopf.v);imag(hopf.v)]),
        parameterbounds=parameterbounds,
        δ=δ,
        δmin=δmin,
        δmax=δmax,
        MaxNumberofSteps = MaxNumberofSteps,
        con_par = nothing,
        NumberOfFails = NumberOfFails
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
