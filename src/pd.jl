struct psol_pd{T}
    profile::Vector{Vector{T}}
    eigenvector::Vector{Vector{T}}
    parameters::Vector{T}
    mesh::Vector{Float64}
    period::T
    ncol::Int
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{ComplexF64,Nothing}
end

function psol_to_pd(pd_guess)
    # construct ns struct
    psol_pd(pd_guess.profile,
        [[0.0]],
        pd_guess.parameters,
        pd_guess.mesh,
        pd_guess.period,
        pd_guess.ncol,
        pd_guess.stability,
        nothing)
end

function pd_q1_approx(jet,pd_guess,τs)
    jac = defsystem_psol_jac(jet,pd_guess,pd_guess,τs; sign1=+)
    q₁ = qr(jac[1:end-1,1:end-1]').Q[:,end]

    # reshape q₁
    γ = pd_guess.profile
    ncol = pd_guess.ncol
    ts = pd_guess.mesh
    ntst = convert(Int,(length(γ)-1)/ncol)
    dims = length(γ[1])
    q₁ = [reshape(q₁[1:dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1]

    # scale q₁ such that ∫ <q₁,q₁> dτ = 1
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    q₁ζ =  interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
    c = sqrt(1/(sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,q₁ζ))/2))
    c*q₁
end 

# for continuation as pd
function pd_w_approx(jet,periodicsolution,τs)
    γ = periodicsolution.profile
    dims = length(γ[1])
    ts = periodicsolution.mesh
    ncol = periodicsolution.ncol
    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

    w = zeros(2dims*(ntst*ncol+1))
    w = reshape(w,4,:)
    w = Vector{eltype(w)}[eachcol(w)...] 
    jac = DDEBifTool.differential_equation_part_complex(jet,w,periodicsolution,τs)
    w = qr([jac[1:end-1,:] rand(2dims*(ntst*ncol+1)); rand(2dims*(ntst*ncol+1)+1)']').Q[:,end]
    w = Vector{eltype(w)}[eachcol(reshape(w[1:end-1],4,:))...] 
    wζ =  interpolate.(colpoints[:],0.0,Ref(w),Ref(ts),ncol)
        
    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    c = sqrt(2.0 / sum(wi .* repeat(weights,ntst) .* dot.(wζ,wζ)))
    @set periodicsolution.eigenvector = c*w
end

function psol_pd_res_as_ns(jet,periodicsolution,psol_ref,τs)
    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    parameters = periodicsolution.parameters
    ω = -π
    w = periodicsolution.eigenvector

    dims = length(γ[1])
    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

     γζ =   interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     wζ =   interpolate.(colpoints[:],0.0,Ref(w),Ref(ts),ncol)
     dγ = d_interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     dw = d_interpolate.(colpoints[:],0.0,Ref(w),Ref(ts),ncol)
    γζs = [hcat(interpolate.(ζ,-τs/T,Ref(γ),Ref(ts),ncol)...) for ζ ∈ colpoints[:]]

    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    dγref = d_interpolate.(colpoints[:],0.0,Ref(psol_ref.profile),T,Ref(ts),ncol)
    phaseconditionI =  sum(wi .* repeat(weights,ntst) .* dot.(γζ, dγref))/2
    partI = [vcat((vec.(dγ)/T - jet.system.(γζs, Ref(parameters))) ...); γ[1] - γ[end]; phaseconditionI ...]

    ϕ₀(τ,θ) =  interpolate(τ,θ,γ,ts,ncol)
    ϕ₁hat(τ) = interpolate(τ,0.0,w,ts,ncol)

    # part II
    @inline R(α) = [cos(α) sin(α); -sin(α) cos(α)]
    res = [zeros(Float64,2dims) for _ ∈ eachindex(ts)]

    M = jet.M
    for i = eachindex(colpoints)
        τ = mod1.(colpoints[i] .- τs/T, 1.0) 
        intervals = searchsortedfirst.(Ref(testintervals), τ) .- 1

        xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γs = [ γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        ws = [ w[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]

        γτ = [L(x,y,τ,ncol) for (x,y,τ) ∈ zip(xs,γs,τ)]
        wτ = [L(x,q,τ,ncol) for (x,q,τ) ∈ zip(xs,ws,τ)]
        dw = dL(xs[1],ws[1],τ[1],ncol)

        res[i] += vec(dw)/T
        for (j,τj) ∈ enumerate(τs[1:end])
           Mj = kron(diagm(ones(dims)), M[j](hcat(γτ...), parameters))
           Rj = kron(R(ω*τj/T), diagm(ones(dims)))
           res[i] -= vec(Mj*Rj*wτ[j])
        end
        res[i] += vec(ω/T*kron([0 -1; 1 0],diagm(ones(dims)))*wτ[1])
    end
    res[end] = vec(w[1] - w[end])
    phaseconditionII =  sum(wi .* repeat(weights,ntst) .* dot.(wζ,wζ))/2 - 1.0

    vcat(partI,vcat(res...),phaseconditionII)
end

function psol_pd_res(jet,q₁,periodicsolution,psol_ref,τs)
    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    parameters = periodicsolution.parameters

    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)

     γζ  =  interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     q₁ζ =  interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
     dγ  = d_interpolate.(colpoints[:],0.0,Ref(γ),Ref(ts),ncol)
     dq₁ = d_interpolate.(colpoints[:],0.0,Ref(q₁),Ref(ts),ncol)
    γζs  = [hcat(interpolate.(ζ,-τs/T,Ref(γ),Ref(ts),ncol)...) for ζ ∈ colpoints[:]]

    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    dγref = d_interpolate.(colpoints[:],0.0,Ref(psol_ref.profile),T,Ref(ts),ncol)
    phaseconditionI =  sum(wi .* repeat(weights,ntst) .* dot.(γζ, dγref))/2
    partI = [vcat((vec.(dγ)/T - jet.system.(γζs, Ref(parameters))) ...); γ[1] - γ[end]; phaseconditionI ...]

    ϕ₀(τ,θ) =  interpolate(τ,θ,γ,ts,ncol)
    ϕ₁hat(τ) = interpolate(τ,0.0,q₁,ts,ncol)

    phaseconditionII =  sum(wi .* repeat(weights,ntst) .* dot.(q₁ζ,q₁ζ))/2 - 1.0
    partII = [vcat(dq₁/T - [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),parameters)*vcat([ϕ₁hat(ζ-τ/T) for τ in τs]...) for ζ in colpoints[:]] ...); q₁[1] + q₁[end]; phaseconditionII ... ]

    vcat(partI,partII)
end

function defsystem_pd_jac(jet,fold_guess,fold_ref,τs)

    M = jet.M

    γ = fold_guess.profile
    # γref = fold_ref.profile
    q₁ = fold_guess.eigenvector
    dims = length(fold_guess.profile[1])
    ncol = fold_guess.ncol
    ntst = convert(Int,(length(γ)-1)/ncol)

    nodes, weights = legendre(ncol)
    ts = fold_guess.mesh
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    par = fold_guess.parameters
    T = fold_guess.period
    testintervals = fold_guess.mesh[1:ncol:end]

    jacI = defsystem_psol_jac(jet,fold_guess,fold_ref,τs,conpar = [1,2])

    jac = zeros(2*dims*(ntst*ncol+1)+3,2*dims*(ntst*ncol+1)+3)
    # derivative partI wrt γ + boundary condition + phase condition
    jac[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)] = jacI[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)]
    # derivative partI wrt parameters and period
    jac[1:dims*(ntst*ncol+1)+1,end-2:end] = jacI[1:dims*(ntst*ncol+1)+1,end-2:end]
    # derivative partII wrt q₁ + boundary condition
    jac[dims*(ntst*ncol+1)+2:end-2,dims*(ntst*ncol+1)+1:end-3] = jacI[1:end-2,1:end-3]
    # change boundary condition for q₁
    jac[end-3:end-2,end-4:end-3] *= -1.0

    ϕ₀(τ,θ) =   interpolate(τ,θ,γ,ts,ncol)
    dϕ₀(τ,θ) = d_interpolate(τ,θ,γ,ts,ncol)
    ϕ₁hat(τ) = interpolate(τ,0.0,q₁,ts,ncol)

    for (i,τ) = enumerate(colpoints)
        ζs = mod.(τ .- τs/T, 1.0) 
        intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1

        if length(unique(intervals)) < 3
            println("Less intervals than three")
            break
        end
        

        xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γζ =   [ γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        q₁ζ = [ q₁[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        u =    [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γζ)]
        q₁ζs = [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,q₁ζ)]
        dq₁ζs = [dL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,q₁ζ)]

        dMj = jet.D2v1(hcat(u ...),par, vcat(-q₁ζs  ...))
        for k=0:ncol
            jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += dMj[:,1:dims]*l0(xs[1],ζs[1],k,ncol)*diagm(ones(dims)) 
        end
        # # delay terms
        for j ∈ 2:length(τs)
            for k=0:ncol
                jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),ncol*dims*(intervals[j]-1) .+ range(k*dims+1,(k+1)*dims)] += dMj[:,(j-1)*dims+1:j*dims]*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims))
            end
        end

        # phase condition 
        ti = xs[1][end] - xs[1][begin]
        for k=0:ncol
            jac[end-1,dims*(ntst*ncol+1) + ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] += ti*weights[mod(i-1,ncol) + 1]*L(xs[1],q₁ζ[1],ζs[1],ncol)*l0(xs[1],ζs[1],k,ncol)
        end

        # derivative with respect to parameters of partII
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims), end-2:end-1] = jet.D11v1(hcat(u ...),par, vcat(-q₁ζs  ...))

        # derivative with respect to period of partII
        # partII = [vcat(dq₁/T - [jet.D1(hcat(ϕ₀.(ζ,-τs/T)...),parameters)*vcat([ϕ₁hat(ζ-τ/T) for τ in τs]...) for ζ in colpoints[:]] ...); 
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),end] = - 1/T^2*dL(xs[1],q₁ζ[1],ζs[1],ncol)
        # delay terms
        for j ∈ 2:length(τs)
            jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*dims+1,i*dims),end] += 1/T^2*τs[j]*dMj[:,(j-1)*dims+1:j*dims]*dL(xs[j],γζ[j],ζs[j],ncol) + M[j](hcat(u...),par)*( - τs[j]/T^2*dq₁ζs[j])
        end
    end

    jac
end

# modified ns_res_jac with ω = pi
function pd_res_jac(jet,periodicsolution,periodicsolution_res,τs)

    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    par = periodicsolution.parameters
    w = periodicsolution.eigenvector
    ω = -π

    dims = length(γ[1])
    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    testintervals = ts[1:ncol:end]

    M = jet.M

    jacI = defsystem_psol_jac(jet,periodicsolution,periodicsolution_res,τs,conpar = [1,2])
    jac = zeros(3*dims*(ntst*ncol+1) + 3, 3*dims*(ntst*ncol+1) + 3)
    # derivative partI wrt γ + boundary condition + phase condition
    jac[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)] = jacI[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)]
    # derivative partI wrt parameters and period
    jac[1:dims*(ntst*ncol+1)+1,end-2:end] = jacI[1:dims*(ntst*ncol+1)+1,end-2:end]

    @inline R(α) = [cos(α) sin(α); -sin(α) cos(α)]
    Rτs = kron(hcat(R.(ω*τs/T)...),diagm(ones(dims)))

    for (i,ζ) = enumerate(colpoints)
        ζs = mod1.(ζ .- τs/T, 1.0) 
        intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1

        xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γs = [ γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        ws = [ w[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γζ =  [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γs)]
        wζs = [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,ws)]
        wζ = wζs[1]
        dwζs = [dL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,ws)]
        dγζs = [dL(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γs)]
        dwζ = dwζs[1]

        for k=0:ncol
            jac[dims*(ntst*ncol+1) + 1 .+ range((i-1)*2*dims+1,i*2*dims),dims*(ntst*ncol+1) .+ ncol*2*dims*(intervals[1]-1) .+ range(k*2*dims+1,(k+1)*2*dims)] += 
                dlj(xs[1],ζs[1],k,ncol)*diagm(ones(2*dims))/T +
                kron(ω/T*[0 -1; 1 0],diagm(ones(dims)))*l0(xs[1],ζs[1],k,ncol)
        end
        # delay terms
        for (j,τ) ∈ enumerate(τs[1:end])
            Mj = kron(diagm(ones(dims)), M[j](hcat(γζ...), par))
            Rj = kron(R(ω*τ/T), diagm(ones(dims)))
            for k=0:ncol
                jac[dims*(ntst*ncol+1) + 1 .+ range((i-1)*2*dims+1,i*2*dims),dims*(ntst*ncol+1) .+ ncol*2*dims*(intervals[j]-1) .+ range(k*2*dims+1,(k+1)*2*dims)] -= 
                    Mj*Rj*l0(xs[j],ζs[j],k,ncol)*diagm(ones(2*dims))
            end
        end

        # derivative with respect to parameters of partII
        # Main.xx[] = wζs
        # return nothing
        v11 = hcat([ cos(ω*τ/T)*wζs[j][1:dims] + sin(ω*τ/T)*wζs[j][dims+1:end] for (j,τ) ∈ enumerate(τs)]...)
        v12 = hcat([-sin(ω*τ/T)*wζs[j][1:dims] + cos(ω*τ/T)*wζs[j][dims+1:end] for (j,τ) ∈ enumerate(τs)]...)
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims), end-2:end-1] = 
            -[jet.D11v1(hcat(γζ ...),par, v11);
              jet.D11v1(hcat(γζ ...),par, v12)]

        # u'/T - f_γ ( τⱼ ->  cos(ω τⱼ/T) u(τ - τⱼ/T) + sin(ω τⱼ/T) v(τ - τⱼ/T) - ω/T v
        # v'/T - f_γ ( τⱼ -> -sin(ω τⱼ/T) u(τ - τⱼ/T) + cos(ω τⱼ/T) v(τ - τⱼ/T) + ω/T u

        dMj_re = jet.D2v1(hcat(γζ...),par,vcat(v11...))
        dMj_im = jet.D2v1(hcat(γζ...),par,vcat(v12...))
        for k=0:ncol
            jac[dims*(ntst*ncol+1) + 1 .+ range((i-1)*2dims+1,i*2dims),ncol*dims*(intervals[1]-1) .+ range(k*dims+1,(k+1)*dims)] -= 
                [dMj_re[:,1:dims]*l0(xs[1],ζs[1],k,ncol)*diagm(ones(dims));
                 dMj_im[:,1:dims]*l0(xs[1],ζs[1],k,ncol)*diagm(ones(dims))]
        end
        # delay terms
        for (j,τ) ∈ enumerate(τs[1:end])
            for k=0:ncol
                jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims),ncol*dims*(intervals[j]-1) .+ range(k*dims+1,(k+1)*dims)] -= 
                [dMj_re[:,(j-1)*dims+1:j*dims]*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims))
                 dMj_im[:,(j-1)*dims+1:j*dims]*l0(xs[j],ζs[j],k,ncol)*diagm(ones(dims))];
            end
        end
        
        # derivative with respect to T of partII
        dv1 = -ω/T^2*vcat([-τ*sin(ω*τ/T)*wζs[j][1:dims] + τ*cos(ω*τ/T)*wζs[j][dims+1:end] for (j,τ) ∈ enumerate(τs)]...)
        dv2 = -ω/T^2*vcat([-τ*cos(ω*τ/T)*wζs[j][1:dims] - τ*sin(ω*τ/T)*wζs[j][dims+1:end] for (j,τ) ∈ enumerate(τs)]...)
        dv1 += 1/T^2*vcat([τ*( cos(ω*τ/T)*dwζs[j][1:dims] + sin(ω*τ/T)*dwζs[j][dims+1:end]) for (j,τ) ∈ enumerate(τs)]...)
        dv2 += 1/T^2*vcat([τ*(-sin(ω*τ/T)*dwζs[j][1:dims] + cos(ω*τ/T)*dwζs[j][dims+1:end]) for (j,τ) ∈ enumerate(τs)]...)
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims), end] = -vec(dwζ)/T^2 + 
            [-jet.D1(hcat(γζ ...),par)*dv1 - -ω*wζ[dims+1:end]/T^2;
             -jet.D1(hcat(γζ ...),par)*dv2 + -ω*wζ[1:dims]/T^2]
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims), end] += 
            [-1/T^2*jet.D2v1(hcat(γζ...),par,vcat(v11...))*vcat([τ*dγζs[j] for (j,τ) ∈ enumerate(τs)]...);
             -1/T^2*jet.D2v1(hcat(γζ...),par,vcat(v12...))*vcat([τ*dγζs[j] for (j,τ) ∈ enumerate(τs)]...)]

        # phaseconditionII = sum(wi .* repeat(weights,ntst) .* dot.(wζ,wζ))/2 - 1.0
        # phaseconditionIII = sum(wi .* repeat(weights,ntst) .* dot.([p[1:dims] for p in wζ],[p[dims+1:end] for p in wζ]))/2
        ti = xs[1][end] - xs[1][begin]
        for k=0:ncol
            jac[end-1,dims*(ntst*ncol+1) .+ ncol*2dims*(intervals[1]-1) .+ range(k*2dims+1,(k+1)*2dims)] += 
                ti*weights[mod(i-1,ncol) + 1]*wζ*l0(xs[1],ζs[1],k,ncol)
        end

    end
    jac[dims*(ntst*ncol+1) + 1 .+ range(2*dims*ntst*ncol+1,2*dims*(ntst*ncol+1)),dims*(ntst*ncol+1) .+ range(1,2*dims)] = diagm(ones(2*dims))
    jac[dims*(ntst*ncol+1) + 1 .+ range(2*dims*ntst*ncol+1,2*dims*(ntst*ncol+1)),end-2*dims+1-3:end-3] = -diagm(ones(2*dims))

    jac
end
