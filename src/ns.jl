struct psol_ns{T}
    profile::Vector{Vector{T}}
    eigenvector::Vector{Vector{T}}
    parameters::Vector{T}
    mesh::Vector{Float64}
    period::T
    ncol::Int
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{ComplexF64,Nothing}
    omega::T
end

function psol_to_ns(ns_guess)
    # find frequency
    indx = findmin(abs.(abs.(ns_guess.stability) .- 1.0))
    ω = atan(imag(ns_guess.stability[indx[2]]),real(ns_guess.stability[indx[2]]))

    # construct ns struct
    psol_ns(ns_guess.profile,
        [[0.0]],
        ns_guess.parameters,
        ns_guess.mesh,
        ns_guess.period,
        ns_guess.ncol,
        ns_guess.stability,
        nothing,
        ω)
end

function psol_ns_res(jet,periodicsolution,psol_ref,τs)
    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    parameters = periodicsolution.parameters
    ω = periodicsolution.omega
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
    phaseconditionIII = sum(wi .* repeat(weights,ntst) .* dot.([p[1:dims] for p in wζ],[p[dims+1:end] for p in wζ]))/2

    vcat(partI,vcat(res...),phaseconditionII,phaseconditionIII)
end

function psol_ns_res_broatcast(jet,periodicsolution,psol_ref,τs)
    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    parameters = periodicsolution.parameters
    ω = periodicsolution.omega
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
    @inline R(α) = [cos(α) -sin(α); sin(α) cos(α)]
    Rτs = kron(hcat(R.(ω*τs/T)...),diagm(ones(dims)))

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
    phaseconditionIII = sum(wi .* repeat(weights,ntst) .* dot.([p[1:dims] for p in wζ],[p[dims+1:end] for p in wζ]))/2

    vcat(partI,vcat(res...),phaseconditionII,phaseconditionIII)
end

function ns_w_approx(jet,periodicsolution,τs)
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
    # w = qr([jac[1:end-1,:] rand(2dims*(ntst*ncol+1)); rand(2dims*(ntst*ncol+1)+1)']').Q[:,end]
    # w = qr(jac[1:end,:]').Q[:,end]
    # @show display([zeros(2dims*(ntst*ncol+1)); 1.0])

    # @show display([jac[1:end-1,:] rand(2dims*(ntst*ncol+1)); rand(2dims*(ntst*ncol+1)+1)'])
    # w = [jac[1:end-1,:] rand(2dims*(ntst*ncol+1)); rand(2dims*(ntst*ncol+1)+1)']\[zeros(2dims*(ntst*ncol+1)); 1.0]
    # w = nullspace(jac[1:end-1,:],atol=1e-02)[:,1]
    # compute nullspace via scd
    _, s, V = svd(jac)
    indxmin = last(findmin(s))
    w = V[:, indxmin]
    @show display(jac[1:end-1,:]*w[1:end])
    w = Vector{eltype(w)}[eachcol(reshape(w[1:end],4,:))...] 
    wζ =  interpolate.(colpoints[:],0.0,Ref(w),Ref(ts),ncol)
        
    testintervals = ts[1:ncol:end]
    wi = repeat(testintervals[2:end] - testintervals[1:end-1],1,ncol)'[:]
    c = sqrt(1.0 / sum(wi .* repeat(weights,ntst) .* dot.(wζ,wζ)))

    @set periodicsolution.eigenvector = c*(w  + [[p[dims+1:end]; -p[1:dims]] for p in w])
end

function differential_equation_part_complex(jet,w,periodicsolution,τs)

    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    par = periodicsolution.parameters
    ω = isdefined(periodicsolution, :omega) ? periodicsolution.omega : -π

    dims = length(γ[1])
    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    testintervals = ts[1:ncol:end]

    M = jet.M

    jac = zeros(2*dims*(ntst*ncol+1) + 1, 2*dims*(ntst*ncol+1))
    @inline R(α) = [cos(α) sin(α); -sin(α) cos(α)]

    for (i,ζ) = enumerate(colpoints)
        ζs = mod1.(ζ .- τs/T, 1.0) 
        intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1

        xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γs = [ γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γζ = [L(x,y,ζ,ncol) for (ζ,x,y) in zip(ζs,xs,γs)]
        wζ = L(xs[1],w[1+(intervals[1]-1)*ncol:1+intervals[1]*ncol],ζs[1],ncol)

        for k=0:ncol
            jac[(i-1)*2*dims+1:i*2*dims,ncol*2*dims*(intervals[1]-1) .+ range(k*2*dims+1,(k+1)*2*dims)] += 
                dlj(xs[1],ζs[1],k,ncol)*diagm(ones(2*dims))/T +
                kron(ω/T*[0 -1; 1 0],diagm(ones(dims)))*l0(xs[1],ζs[1],k,ncol)
        end
        # delay terms
        for (j,τ) ∈ enumerate(τs[1:end])
            Mj = kron(diagm(ones(dims)), M[j](hcat(γζ...), par))
            Rj = kron(R(ω*τ/T), diagm(ones(dims)))
            for k=0:ncol
                jac[(i-1)*2*dims+1:i*2*dims,ncol*2*dims*(intervals[j]-1) .+ range(k*2*dims+1,(k+1)*2*dims)] -= 
                    Mj*Rj*l0(xs[j],ζs[j],k,ncol)*diagm(ones(2*dims))
            end
        end

        # phase condition I
        ti = xs[1][end] - xs[1][begin]
        for k=0:ncol
            jac[end,ncol*2dims*(intervals[1]-1) .+ range(k*2dims+1,(k+1)*2dims)] += 
                ti*weights[mod(i-1,ncol) + 1]*wζ*l0(xs[1],ζs[1],k,ncol)
        end
    end
    jac[2*dims*ntst*ncol+1:2*dims*(ntst*ncol+1),1:2*dims] = diagm(ones(2*dims))
    jac[2*dims*ntst*ncol+1:2*dims*(ntst*ncol+1),end-2*dims+1:end] = -diagm(ones(2*dims))

    jac
end

function ns_res_jac(jet,periodicsolution,periodicsolution_res,τs)

    γ = periodicsolution.profile
    T = periodicsolution.period
    ncol = periodicsolution.ncol
    ts = periodicsolution.mesh
    par = periodicsolution.parameters
    ω = periodicsolution.omega
    w = periodicsolution.eigenvector

    dims = length(γ[1])
    ntst = convert(Int,(length(γ)-1)/ncol)
    nodes,weights = legendre(ncol)
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    testintervals = ts[1:ncol:end]

    M = jet.M

    jacI = defsystem_psol_jac(jet,periodicsolution,periodicsolution_res,τs,conpar = [1,2])
    jac = zeros(3*dims*(ntst*ncol+1) + 4, 3*dims*(ntst*ncol+1) + 4)
    # derivative partI wrt γ + boundary condition + phase condition
    jac[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)] = jacI[1:dims*(ntst*ncol+1)+1,1:dims*(ntst*ncol+1)]
    # derivative partI wrt parameters and period
    jac[1:dims*(ntst*ncol+1)+1,end-3:end-1] = jacI[1:dims*(ntst*ncol+1)+1,end-2:end]

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
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims), end-3:end-2] = 
            -[jet.D11v1(hcat(γζ ...),par, v11);
              jet.D11v1(hcat(γζ ...),par, v12)]

        # u'/T - f_γ ( τⱼ ->  cos(ω τⱼ/T) u(τ - τⱼ/T) + sin(ω τⱼ/T) v(τ - τⱼ/T) - ω/T v == 0
        # v'/T - f_γ ( τⱼ -> -sin(ω τⱼ/T) u(τ - τⱼ/T) + cos(ω τⱼ/T) v(τ - τⱼ/T) + ω/T u == 0
        #
        # note that we have two degrees of freedom in (u,v), namely
        # 1.  (u,v) -> c_1 (u,v)
        # 2.  (u,v) -> c_2 (v,-u)
        #
        # v'/T - f_γ ( τⱼ ->  cos(ω τⱼ/T) v(τ - τⱼ/T) - sin(ω τⱼ/T) u(τ - τⱼ/T) + ω/T u == 0
        # u'/T - f_γ ( τⱼ ->  sin(ω τⱼ/T) v(τ - τⱼ/T) + cos(ω τⱼ/T) u(τ - τⱼ/T) - ω/T v == 0
        #
        # thus if (u,v) is a solutions, then so is c1 (u,v) + c2 (v,-u)
        # 
        # We can use this freedom two satisfy the following normalization conditions
        # 1. ∫ <u,u> + <v,v>  dτ - 1 == 0
        # 2. ∫ <u,v> = 0
        # 
        # Indeed, we have that
        #
        # 1. (c1^2 + c2^2) ∫ <u,u> + <v,v> dτ - 1 == 0
        # 2. (c1^2 - c2^2) ∫ <u,v> dτ  = 0
        #
        # Thus let c1^2 == c2^2, so that the second intergral vanishes. This yields that
        # 1. 2 c1^2  ∫ <u,u> + <v,v> dτ - 1 == 0
        # 
        # Solveing for c1 yields
        # 1. c1  == sqrt(1/∫ <u,u> + <v,v> dτ/2) (Nothing that c1 > 0)

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
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims), end-1] = -vec(dwζ)/T^2 + 
            [-jet.D1(hcat(γζ ...),par)*dv1 - -ω*wζ[dims+1:end]/T^2;
             -jet.D1(hcat(γζ ...),par)*dv2 + -ω*wζ[1:dims]/T^2]
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims), end-1] += 
            [-1/T^2*jet.D2v1(hcat(γζ...),par,vcat(v11...))*vcat([τ*dγζs[j] for (j,τ) ∈ enumerate(τs)]...);
             -1/T^2*jet.D2v1(hcat(γζ...),par,vcat(v12...))*vcat([τ*dγζs[j] for (j,τ) ∈ enumerate(τs)]...)]

        # derivative with respect to ω of partII
        dv1 = 1/T*vcat([-τ*sin(ω*τ/T)*wζs[j][1:dims] + τ*cos(ω*τ/T)*wζs[j][dims+1:end] for (j,τ) ∈ enumerate(τs)]...)
        dv2 = 1/T*vcat([-τ*cos(ω*τ/T)*wζs[j][1:dims] - τ*sin(ω*τ/T)*wζs[j][dims+1:end] for (j,τ) ∈ enumerate(τs)]...)
        jac[dims*(ntst*ncol+1)+1 .+ range((i-1)*2dims+1,i*2dims), end] = 
            [-jet.D1(hcat(γζ ...),par)*dv1 - wζ[dims+1:end]/T;
             -jet.D1(hcat(γζ ...),par)*dv2 + wζ[1:dims]/T]

        # phaseconditionII = sum(wi .* repeat(weights,ntst) .* dot.(wζ,wζ))/2 - 1.0
        # phaseconditionIII = sum(wi .* repeat(weights,ntst) .* dot.([p[1:dims] for p in wζ],[p[dims+1:end] for p in wζ]))/2
        ti = xs[1][end] - xs[1][begin]
        for k=0:ncol
            jac[end-2,dims*(ntst*ncol+1) .+ ncol*2dims*(intervals[1]-1) .+ range(k*2dims+1,(k+1)*2dims)] += 
                ti*weights[mod(i-1,ncol) + 1]*wζ*l0(xs[1],ζs[1],k,ncol)
            jac[end-1,dims*(ntst*ncol+1) .+ ncol*2dims*(intervals[1]-1) .+ range(k*2dims+1,k*2dims+dims)] += 
                ti*weights[mod(i-1,ncol) + 1]*wζ[dims+1:end]*l0(xs[1],ζs[1],k,ncol)/2
            jac[end-1,dims*(ntst*ncol+1) .+ ncol*2dims*(intervals[1]-1) .+ range(k*2dims+dims+1,(k+1)*2dims)] += 
                ti*weights[mod(i-1,ncol) + 1]*wζ[1:dims]*l0(xs[1],ζs[1],k,ncol)/2
        end

    end
    jac[dims*(ntst*ncol+1) + 1 .+ range(2*dims*ntst*ncol+1,2*dims*(ntst*ncol+1)),dims*(ntst*ncol+1) .+ range(1,2*dims)] = diagm(ones(2*dims))
    jac[dims*(ntst*ncol+1) + 1 .+ range(2*dims*ntst*ncol+1,2*dims*(ntst*ncol+1)),end-2*dims+1-4:end-4] = -diagm(ones(2*dims))

    jac
end
