struct psol{T}
    profile::Vector{Vector{T}}
    parameters::Vector{T}
    mesh::Vector{Float64}
    period::T
    ncol::Int
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{ComplexF64,Nothing}
end

struct psol_fold{T}
    profile::Vector{Vector{T}}
    eigenvector::Vector{Vector{T}}
    parameters::Vector{T}
    mesh::Vector{Float64}
    period::T
    ncol::Int
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{ComplexF64,Nothing}
    beta::T
end

function psol_res(jet,periodicsolution,psol_ref,τs)
    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    parameters = periodicsolution.parameters

    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

    # dγ = [hcat(d_interpolate(ζ,0.0,γ,T,ts,ncol)...) for ζ ∈ colpoints[:]]
     γζ =   interpolate.(colpoints[:],0.0,Ref(γ),T,Ref(ts),ncol)
     dγ = d_interpolate.(colpoints[:],0.0,Ref(γ),T,Ref(ts),ncol)
    γζs = [hcat(interpolate.(ζ,-τs,Ref(γ),T,Ref(ts),ncol)...) for ζ ∈ colpoints[:]]

    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    dγref = d_interpolate.(colpoints[:],0.0,Ref(psol_ref.profile),T,Ref(ts),ncol)
    phasecondition =  sum(wi .* repeat(weights,ntst) .* dot.(γζ, dγref))/2
    [vcat((vec.(dγ)/T - jet.system.(γζs, Ref(parameters))) ...); γ[1] - γ[end]; phasecondition ...]
end

function +(p1::psol,p2::psol)
    psol(
         p1.profile .+ p2.profile,
         p1.parameters + p2.parameters,
         p2.mesh,
         p1.period + p2.period,
         p2.ncol,
         nothing,
         nothing
    )
end

function -(p1::psol,p2::psol)
    psol(
         p1.profile .-p2.profile,
         p1.parameters - p2.parameters,
         p2.mesh,
         p1.period - p2.period,
         p2.ncol,
         nothing,
         nothing
    )
end

function *(δ,p::psol)
    psol(
         δ*p.profile,
         δ*p.parameters,
         p.mesh,
         δ*p.period,
         p.ncol,
         nothing,
         nothing
    )
end 

function +(p1::psol_fold,p2::psol_fold)
    psol_fold(
         p1.profile .+ p2.profile,
         p1.eigenvector .+ p2.eigenvector,
         p1.parameters + p2.parameters,
         p2.mesh,
         p1.period + p2.period,
         p2.ncol,
         nothing,
         nothing,
         p1.beta + p2.beta
    )
end

function -(p1::psol_fold,p2::psol_fold)
    psol_fold(
         p1.profile .-p2.profile,
         p1.eigenvector .- p2.eigenvector,
         p1.parameters - p2.parameters,
         p2.mesh,
         p1.period - p2.period,
         p2.ncol,
         nothing,
         nothing,
         p1.beta - p2.beta
    )
end

function *(δ,p::psol_fold)
    psol_fold(
         δ*p.profile,
         δ*p.eigenvector,
         δ*p.parameters,
         p.mesh,
         δ*p.period,
         p.ncol,
         nothing,
         nothing,
         δ*p.beta
    )
end 

function remesh_psol(psol1, new_mesh, new_ncol)
    γ = psol1.profile
    ncol = psol1.ncol
    ts = psol1.mesh 
    γ_new = [[γ[1]]; interpolate.(new_mesh[2:end-1],0.0,Ref(γ),Ref(ts),ncol); [γ[end]]]
    psol(γ_new, psol1.parameters, collect(new_mesh), psol1.period, new_ncol, nothing, nothing)
end

function fold_q1_approx(jet,fold_guess,τs)
    γ = fold_guess.profile
    dims = length(γ[1])
    ncol = fold_guess.ncol
    ntst = convert(Int,(length(γ)-1)/ncol)

    T = fold_guess.period
    ts = fold_guess.mesh
    testintervals = fold_guess.mesh[1:ncol:end]
    par = fold_guess.parameters
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

    ϕ₀(τ,θ) =   interpolate(τ,θ,γ,ts,ncol)
    dγ(τ,θ) = d_interpolate(τ,θ,γ,ts,ncol)

    RHS = dγ.(colpoints[:],0.0)/T + [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),par)*vcat([τ*dγ(ζ,-τ/T) for τ in τs]...) for ζ in colpoints[:]]/T
    Bex = [vcat(RHS...); zeros(dims+1)]
    jac = defsystem_psol_jac(jet,fold_guess,fold_guess,τs)
    q1,_ = borderedInverse(jac[1:end-1,1:end-2],-Bex,ntst,ncol,dims,ts,normalization=true);
    q1
end 

function fold_tangent(jet,fold_guess,τs)
    # jac = defsystem_psol_jac(jet,fold_guess,fold_guess,τs)
    # v = qr(jac[1:end,1:end .!= end-1]').Q[:,end-1]
    q₁ = fold_q1_approx(jet,fold_guess,τs)
    β = 1.0
    v  = [q₁...; β]

    γ = fold_guess.profile
    ncol = fold_guess.ncol
    ts = fold_guess.mesh

    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

    # q₁ = [reshape(v[1:dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1]
    # β = v[end]

    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    q₁ζ =   interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
    vβnorm =  sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,q₁ζ))/2 + β^2
    v /= sqrt(vβnorm)
end 

function borderedInverse(jac,rhs,ntst,ncol,dims,ts;normalization=true)
    jacex = [jac[1:end-1,:] rand(typeof(first(jac)),dims*(ntst*ncol+1))
             rand(typeof(first(jac)),1,dims*(ntst*ncol+1)) 0.0]
    p1 = jacex'\[zeros(dims*(ntst*ncol+1)); 1.0]
    if normalization
        sol = [jac p1]\rhs
    else
        q1 = jacex\[zeros(dims*(ntst*ncol+1)); 1.0]
        sol = [jac[1:end-1,:] p1[1:end-1]; q1[1:end-1]' 0.0]\rhs
        # sol = jacex\rhs
    end
    sol = [vec(reshape(sol[dims*(i-1)+1:dims*i],dims,1)) for i ∈ eachindex(ts)]
    sol, p1
end

function psol_fold_res(jet,q₁,periodicsolution,psol_ref,τs)
    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    parameters = periodicsolution.parameters

    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

    # dγ = [hcat(d_interpolate(ζ,0.0,γ,T,ts,ncol)...) for ζ ∈ colpoints[:]]
     γζ =   interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     q₁ζ =   interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
     dγ = d_interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     dq₁ = d_interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
    γζs = [hcat(interpolate.(ζ,-τs/T,Ref(γ),Ref(ts),ncol)...) for ζ ∈ colpoints[:]]

    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    dγref = d_interpolate.(colpoints[:],0.0,Ref(psol_ref.profile),T,Ref(ts),ncol)
    phaseconditionI =  sum(wi .* repeat(weights,ntst) .* dot.(γζ, dγref))/2
    # phasecondition = sum(dot.(γζ,dγ))
    partI = [vcat((vec.(dγ)/T - jet.system.(γζs, Ref(parameters))) ...); γ[1] - γ[end]; phaseconditionI ...]

    ϕ₀(τ,θ) =   interpolate(τ,θ,γ,ts,ncol)
    dϕ₀(τ,θ) = d_interpolate(τ,θ,γ,ts,ncol)
    ϕ₁hat(τ) = interpolate(τ,0.0,q₁,ts,ncol)

    # dγ/T + dq₁/T + [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),par)*vcat([τ*dγ(ζ,-τ/T) - ϕ₁hat(ζ-τ/T) for τ in τs]...) for ζ in colpoints[:]]/T

    phaseconditionII =  sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,dγref))/2
    partII = [vcat(dγ/T + dq₁/T + [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),parameters)*vcat([τ*dϕ₀(ζ,-τ/T)/T - ϕ₁hat(ζ-τ/T) for τ in τs]...) for ζ in colpoints[:]] ...); q₁[1] - q₁[end]; phaseconditionII ... ]

    vcat(partI,partII)
end

# function psol_pd_res(jet,q₁,periodicsolution,psol_ref,τs)
#     γ = periodicsolution.profile
#     T = periodicsolution.period
#     ncol = periodicsolution.ncol
#     ts = periodicsolution.mesh
#     parameters = periodicsolution.parameters
# 
#     ntst = convert(Int,(length(γ)-1)/ncol)
#     nodes,weights = legendre(ncol)
#     colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
# 
#      γζ  =  interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
#      q₁ζ =  interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
#      dγ  = d_interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
#      dq₁ = d_interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
#     γζs  = [hcat(interpolate.(ζ,-τs/T,Ref(γ),Ref(ts),ncol)...) for ζ ∈ colpoints[:]]
# 
#     testintervals = ts[1:ncol:end]
#     wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
#     dγref = d_interpolate.(colpoints[:],0.0,Ref(psol_ref.profile),T,Ref(ts),ncol)
#     phaseconditionI =  sum(wi .* repeat(weights,ntst) .* dot.(γζ, dγref))/2
#     partI = [vcat((vec.(dγ)/T - jet.system.(γζs, Ref(parameters))) ...); γ[1] - γ[end]; phaseconditionI ...]
# 
#     ϕ₀(τ,θ) =  interpolate(τ,θ,γ,ts,ncol)
#     ϕ₁hat(τ) = interpolate(τ,0.0,q₁,ts,ncol)
# 
#     phaseconditionII =  sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,q₁ζ))/2 - 1.0
#     partII = [vcat(dq₁/T - [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),parameters)*vcat([ϕ₁hat(ζ-τ/T) for τ in τs]...) for ζ in colpoints[:]] ...); q₁[1] + q₁[end]; phaseconditionII ... ]
# 
#     vcat(partI,partII)
# end

function psol_fold_res_β(jet,q₁,β,periodicsolution,psol_ref,τs)
    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    parameters = periodicsolution.parameters

    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

     γζ  =   interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     q₁ζ =   interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
     dγ = d_interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     dq₁ = d_interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
    γζs = [hcat(interpolate.(ζ,-τs/T,Ref(γ),Ref(ts),ncol)...) for ζ ∈ colpoints[:]]

    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    dγref = d_interpolate.(colpoints[:],0.0,Ref(psol_ref.profile),T,Ref(ts),ncol)
    phaseconditionI =  sum(wi .* repeat(weights,ntst) .* dot.(γζ, dγref))/2
    # phasecondition = sum(dot.(γζ,dγ))
    partI = [vcat((vec.(dγ)/T - jet.system.(γζs, Ref(parameters))) ...); γ[1] - γ[end]; phaseconditionI ...]

    ϕ₀(τ,θ) =   interpolate(τ,θ,γ,ts,ncol)
    dϕ₀(τ,θ) = d_interpolate(τ,θ,γ,ts,ncol)
    ϕ₁hat(τ) = interpolate(τ,0.0,q₁,ts,ncol)

    phaseconditionII =  sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,dγref))/2
    # phaseconditionII =  sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,dq₁))/2
    phaseconditionIII =  sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,q₁ζ))/2 + β^2 - 1.0
    partII = [vcat(β*dγ/T + dq₁/T + [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),parameters)*vcat([β*τ*dϕ₀(ζ,-τ/T)/T - ϕ₁hat(ζ-τ/T) for τ in τs]...) for ζ in colpoints[:]] ...); 
        q₁[1] - q₁[end]; phaseconditionII; phaseconditionIII ... ]
    vcat(partI,partII)
end

function defsystem_psol_jac(jet,psol_guess,psol_ref,τs; conpar = 2, sign1=-)
    M = jet.M

    γ = psol_guess.profile
    γref = psol_ref.profile
    dims = length(psol_guess.profile[1])
    ncol = psol_guess.ncol
    ntst = convert(Int,(length(γ)-1)/ncol)

    nodes, weights = legendre(ncol)
    ts = psol_guess.mesh
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    par = psol_guess.parameters
    T = psol_guess.period
    testintervals = psol_guess.mesh[1:ncol:end]

    Jac = zeros(dims*(ntst*ncol+1) + 2,dims*(ntst*ncol+1) + 1 + length(conpar))
    for (i,τ) = enumerate(colpoints)
        ζs = mod1.(τ .- τs/T, 1.0) 
        intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1

        xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γζ = [ γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γrefζ = [ γref[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        u = [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]

        for k=0:ncol
            Jac[(i-1)*dims+1:i*dims,ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += 1/T*dlj(xs[1],ζs[1],k,ncol)*diagm(ones(dims)) - M[1](hcat(u...), par)*l0(xs[1],ζs[1],k,ncol)*diagm(ones(dims)) 
        end
        # delay terms
        for j ∈ 2:length(τs)
            for k=0:ncol
                Jac[(i-1)*dims+1:i*dims,ncol*dims*(intervals[j]-1) .+ range(k*dims+1,(k+1)*dims)] -= M[j](hcat(u...),par)*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims)) 
            end
        end
        # phase condition
        ti = xs[1][end] - xs[1][begin]
        for k=0:ncol
            Jac[end-1,ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += 0.5*ti*weights[mod(i-1,ncol) + 1]*dL(xs[1],γrefζ[1],ζs[1],ncol)*l0(xs[1],ζs[1],k,ncol)
        end

        # derivative with respect to parameter
        Jac[(i-1)*dims+1:i*dims,dims*(ntst*ncol+1)+1:dims*(ntst*ncol+1)+length(conpar)] = -jet.D01(hcat(u...), par)[:,conpar]

        # derivative with respect to period
        Jac[(i-1)*dims+1:i*dims,end] = -1/T^2*dL(xs[1],γζ[1],ζs[1],ncol)
        # delay terms
        for j ∈ 2:length(τs)
            Jac[(i-1)*dims+1:i*dims,end] -= 1/T^2*τs[j]*M[j](hcat(u...),par)*dL(xs[j],γζ[j],ζs[j],ncol)
        end

    end
    Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1),1:dims] = diagm(ones(dims))
    Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1),end-dims-(1+length(conpar))+1:end-(1+length(conpar))] = sign1(diagm(ones(dims)))

    Jac
end

function defsystem_fold_jac(jet,fold_guess,fold_ref,τs)

    M = jet.M

    γ = fold_guess.profile
    # γref = fold_ref.profile
    q₁ = fold_guess.eigenvector
    β  = fold_guess.beta
    dims = length(fold_guess.profile[1])
    ncol = fold_guess.ncol
    ntst = convert(Int,(length(γ)-1)/ncol)

    nodes, weights = legendre(ncol)
    ts = fold_guess.mesh
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    par = fold_guess.parameters
    T = fold_guess.period
    testintervals = fold_guess.mesh[1:ncol:end]

    dγ = d_interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)

    jacI = defsystem_psol_jac(jet,fold_guess,fold_guess,τs,conpar = [1,2])

    jac_sym = zeros(2*dims*(ntst*ncol+1)+4,2*dims*(ntst*ncol+1)+4)
    # derivative partI wrt γ + boundary condition + phase condition
    jac_sym[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)] = jacI[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)]
    # derivative partI wrt parameters and period
    jac_sym[1:dims*(ntst*ncol+1)+1,end-3:end-1] = jacI[1:dims*(ntst*ncol+1)+1,end-2:end]
    # derivative partII wrt q₁ + boundary condition + phaseconditionII
    jac_sym[dims*(ntst*ncol+1)+2:end-2,dims*(ntst*ncol+1)+1:end-4] = jacI[1:end-1,1:end-3]

    ϕ₀(τ,θ) =   interpolate(τ,θ,γ,ts,ncol)
    dϕ₀(τ,θ) = d_interpolate(τ,θ,γ,ts,ncol)
    ϕ₁hat(τ) = interpolate(τ,0.0,q₁,ts,ncol)

    # derivative partII wrt β
    jac_sym[dims*(ntst*ncol+1)+2:end-dims-3,end] = vcat(dγ/T + [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),par)*vcat([τ*dϕ₀(ζ,-τ/T)/T for τ in τs]...) for ζ in colpoints[:]] ...)
    jac_sym[end-1,end] = 2*fold_guess.beta

    for (i,τ) = enumerate(colpoints)
        ζs = mod.(τ .- τs/T, 1.0) 
        intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1

        xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γζ =   [ γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        q₁ζ = [ q₁[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        u =    [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]
        du =  [dL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]
        ddu =  [ddL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]
        q₁ζs = [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,q₁ζ)]
        dq₁ζs = [dL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,q₁ζ)]

        dMj = jet.D2v1(hcat(u ...),par, vcat(β/T*τs.*du - q₁ζs  ...))
        # dMj = jet.D2v1(hcat(ϕ₀.(τ,-τs/T)...),par, vcat([β*τj*dϕ₀(τ,-τj/T)/T - ϕ₁hat(τ-τj/T) for τj in τs] ...))
        # dMj = zeros(dims,dims*length(τs))

        for k=0:ncol
            jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += β/T*dlj(xs[1],ζs[1],k,ncol)*diagm(ones(dims)) + dMj[:,1:dims]*l0(xs[1],ζs[1],k,ncol)*diagm(ones(dims)) 
        end
        # # delay terms
        for j ∈ 2:length(τs)
            for k=0:ncol
                jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[j]-1) .+ range(k*dims+1,(k+1)*dims)] += β/T*τs[j]*M[j](hcat(u...),par)*dlj(xs[j],ζs[j],k,ncol)*diagm(ones(dims)) + dMj[:,(j-1)*dims+1:j*dims]*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims))
                # jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[j]-1) .+ range(k*dims+1,(k+1)*dims)] += dMj[:,(j-1)*dims+1:j*dims]*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims))
            end
        end

        # phase condition III
        ti = xs[1][end] - xs[1][begin]
        for k=0:ncol
            jac_sym[end-1,dims*(ntst*ncol+1) + ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += ti*weights[mod(i-1,ncol) + 1]*L(xs[1],q₁ζ[1],ζs[1],ncol)*l0(xs[1],ζs[1],k,ncol)
        end

        # derivative with respect to parameters of partII
        jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims), end-3:end-2] = jet.D11v1(hcat(u ...),par, vcat(β/T*τs.*du - q₁ζs  ...))

        # derivative with respect to period of partII
        jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),end-1] = -β/T^2*dL(xs[1],γζ[1],ζs[1],ncol) - 1/T^2*dL(xs[1],q₁ζ[1],ζs[1],ncol)
        # delay terms
        for j ∈ 2:length(τs)
            jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),end-1] += 1/T^2*τs[j]*dMj[:,(j-1)*dims+1:j*dims]*dL(xs[j],γζ[j],ζs[j],ncol) + M[j](hcat(u...),par)*( β*τs[j]*((-1/T^2)*du[j] + τs[j]/T^2*ddu[j]/T) - τs[j]/T^2*dq₁ζs[j])
        end

    end

    jac_sym
end

function defsystem_fold_jac_β(jet,fold_guess,fold_ref,τs)

    M = jet.M

    γ = fold_guess.profile
    # γref = fold_ref.profile
    q₁ = fold_guess.eigenvector
    β  = fold_guess.beta
    dims = length(fold_guess.profile[1])
    ncol = fold_guess.ncol
    ntst = convert(Int,(length(γ)-1)/ncol)

    nodes, weights = legendre(ncol)
    ts = fold_guess.mesh
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    par = fold_guess.parameters
    T = fold_guess.period
    testintervals = fold_guess.mesh[1:ncol:end]

    dγ = d_interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)

    jacI = defsystem_psol_jac(jet,fold_guess,fold_guess,τs,conpar = [1,2])

    jac_sym = zeros(2*dims*(ntst*ncol+1)+4,2*dims*(ntst*ncol+1)+4)
    # derivative partI wrt γ + boundary condition + phase condition
    jac_sym[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)] = jacI[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)]
    # derivative partI wrt parameters and period
    jac_sym[1:dims*(ntst*ncol+1)+1,end-3:end-1] = jacI[1:dims*(ntst*ncol+1)+1,end-2:end]
    # derivative partII wrt q₁ + boundary condition + phaseconditionII
    jac_sym[dims*(ntst*ncol+1)+2:end-2,dims*(ntst*ncol+1)+1:end-4] = jacI[1:end-1,1:end-3]

    ϕ₀(τ,θ) =   interpolate(τ,θ,γ,ts,ncol)
    dϕ₀(τ,θ) = d_interpolate(τ,θ,γ,ts,ncol)
    ϕ₁hat(τ) = interpolate(τ,0.0,q₁,ts,ncol)

    # derivative partII wrt β
    jac_sym[dims*(ntst*ncol+1)+2:end-dims-3,end] = vcat(dγ/T + [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),par)*vcat([τ*dϕ₀(ζ,-τ/T)/T for τ in τs]...) for ζ in colpoints[:]] ...)
    jac_sym[end-1,end] = 2*fold_guess.beta

    for (i,τ) = enumerate(colpoints)
        ζs = mod.(τ .- τs/T, 1.0) 
        intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1

        xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γζ =   [ γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        q₁ζ = [ q₁[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        u =    [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]
        du =  [dL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]
        ddu =  [ddL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]
        q₁ζs = [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,q₁ζ)]
        dq₁ζs = [dL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,q₁ζ)]

        dMj = jet.D2v1(hcat(u ...),par, vcat(β/T*τs.*du - q₁ζs  ...))
        # dMj = jet.D2v1(hcat(ϕ₀.(τ,-τs/T)...),par, vcat([β*τj*dϕ₀(τ,-τj/T)/T - ϕ₁hat(τ-τj/T) for τj in τs] ...))
        # dMj = zeros(dims,dims*length(τs))

        for k=0:ncol
            jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += β/T*dlj(xs[1],ζs[1],k,ncol)*diagm(ones(dims)) + dMj[:,1:dims]*l0(xs[1],ζs[1],k,ncol)*diagm(ones(dims)) 
        end
        # # delay terms
        for j ∈ 2:length(τs)
            for k=0:ncol
                jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[j]-1) .+ range(k*dims+1,(k+1)*dims)] += β/T*τs[j]*M[j](hcat(u...),par)*dlj(xs[j],ζs[j],k,ncol)*diagm(ones(dims)) + dMj[:,(j-1)*dims+1:j*dims]*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims))
                # jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[j]-1) .+ range(k*dims+1,(k+1)*dims)] += dMj[:,(j-1)*dims+1:j*dims]*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims))
            end
        end

        # phase condition III
        ti = xs[1][end] - xs[1][begin]
        for k=0:ncol
            jac_sym[end-1,dims*(ntst*ncol+1) + ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += ti*weights[mod(i-1,ncol) + 1]*L(xs[1],q₁ζ[1],ζs[1],ncol)*l0(xs[1],ζs[1],k,ncol)
        end

        # derivative with respect to parameters of partII
        jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims), end-3:end-2] = jet.D11v1(hcat(u ...),par, vcat(β/T*τs.*du - q₁ζs  ...))

        # derivative with respect to period of partII
        jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),end-1] = -β/T^2*dL(xs[1],γζ[1],ζs[1],ncol) - 1/T^2*dL(xs[1],q₁ζ[1],ζs[1],ncol)
        # delay terms
        for j ∈ 2:length(τs)
            jac_sym[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),end-1] += 1/T^2*τs[j]*dMj[:,(j-1)*dims+1:j*dims]*dL(xs[j],γζ[j],ζs[j],ncol) + M[j](hcat(u...),par)*( β*τs[j]*((-1/T^2)*du[j] + τs[j]/T^2*ddu[j]/T) - τs[j]/T^2*dq₁ζs[j])
        end

    end

    jac_sym
end
