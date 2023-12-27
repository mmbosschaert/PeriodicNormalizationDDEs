struct psol{T}
    profile::Vector{Vector{T}}
    parameters::Vector{T}
    mesh::Vector{Float64}
    period::T
    ncol::Int
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{ComplexF64,Nothing}
end

struct psol_pd{T}
    profile::Vector{Vector{T}}
    eigenvector::Vector{Vector{T}}
    parameters::Vector{T}
    mesh::Vector{Float64}
    period::T
    ncol::Int
    stability::Union{Vector{ComplexF64},Nothing}
    nmfm::Union{Float64,Nothing}
    omega::T
end


function vec(p::psol,con_par)
    [vcat(p.profile...); p.parameters[con_par]; p.period ...]
end

function vec_to_point(vec,point_prev::psol,con_par)
    psol(
            [reshape(vec[1:end-2],2,:)[:,j] for j in 1:length(point_prev.mesh)], 
            [point_prev.parameters[1:con_par-1]; vec[end-1]; point_prev.parameters[con_par+1:end]], 
             point_prev.mesh, vec[end], point_prev.ncol, nothing, nothing)
end

function SetupPsolBranch(jet, con_par, hopf_point::Hopf, τs; parameterbounds=nothing, δ=.001, δmin=1e-06, δmax=0.01, MaxNumberofSteps = 250, ntst = 20, ncol = 3, NumberOfFails = 4)

    # define fine mesh
    t = range(0.0,1.0, ntst*ncol + 1)
    profile = 1e-03*imag([hopf_point.v*s for s in exp.(im*2π*t)])
    T = 2pi/abs(hopf_point.ω)
    psol_guess = psol(profile, hopf_point.parameters, collect(t), T, ncol, nothing, nothing)
    psol_res(jet,psol_guess,psol_guess,τs)

    f = (x,psol1) -> vcat([
        psol_res(jet,
            psol([reshape(x[1:end-2],2,:)[:,j] for j in 1:length(psol1.mesh)],
            [psol1.parameters[1:con_par-1]; x[end-1]; psol1.parameters[con_par+1:end]], psol1.mesh, x[end], 
                psol1.ncol, nothing, nothing),psol1,τs)...; 0.0] ...)

    df = (x,psol1) -> defsystem_psol_jac(jet,
            psol([reshape(x[1:end-2],2,:)[:,j] for j in 1:length(psol1.mesh)],
            [psol1.parameters[1:con_par-1]; x[end-1]; psol1.parameters[con_par+1:end]], psol1.mesh, x[end], 
                psol1.ncol, nothing, nothing),psol1,τs)

    x₀ = vec(psol_guess,con_par)
    jac = df(x₀,psol_guess)
    V = [jac[1:end-1,:]; rand(length(x₀))']\[zeros(length(x₀)-1); 1.0]

    x₀new, _, V_new = newton(f, df, x₀, V, psol_guess; tol=1e-8)
    psol_corrected = vec_to_point(x₀new,psol_guess,con_par)

    # continuation with own newton
    psol_branch = (points = psol[],
        tangents = [],
        stepsizes = [],
        f = f,
        df = df,
        parameterbounds = parameterbounds,
        δ=δ,
        δmin=δmin,
        δmax=δmax,
        MaxNumberofSteps = MaxNumberofSteps,
        con_par = con_par,
        NumberOfFails = NumberOfFails)
    push!(psol_branch.points, psol_corrected)
    push!(psol_branch.tangents, V_new)
    push!(psol_branch.stepsizes, 0.0)
    psol_branch
end

# branch off near psol to psol with double period
function SetupPsolBranch(jet, con_par, psol1::psol_pd, τs; parameterbounds=nothing, δ=.001, δmin=1e-06, δmax=0.01, MaxNumberofSteps = 250, NumberOfFails = 4)

    dims = length(psol1.profile[1])
    ncol = psol1.ncol
    ntst = convert(Int,(length(psol1.profile)-1)/ncol)

    # define fine mesh
    # compute eigenvector
    jac = defsystem_psol_jac_no_mod(jet,psol1,τs)
    s1,s2=size(jac)
    n_ext=s2-s1
    n2_ext=max(n_ext-s1,0)
    B=blockdiag(sparse(1.0I, n2_ext, n2_ext),sparse(jac[1:s1,n_ext+1:s2]))
    Atop=[spzeros(n2_ext,s1) sparse(1.0I, n2_ext, n2_ext)]
    A=[Atop[:,1:n_ext];-jac[1:s1,1:n_ext]]
    Ptau=[spzeros(n_ext,max(0,2*s1-s2)) sparse(1.0I, n_ext, n_ext)]

    M = B\Matrix(A*sparse(1.0I, n_ext, n_ext))
    vals,vecs = eigen(Ptau*M)
    indx = last(findmin(abs.(vals) .+ 1))

    MPfun(x) = vcat(x,(B[end-s1+1:end,:]\(A*x)))
    x = real(MPfun(vecs[:,indx]))
    [c[:] for c in eachcol(reshape(x[1:dims*(ntst*ncol+1)],dims,:))]

    deviation=[c[:] for c in eachcol(reshape(x[1:dims*(ntst*ncol+1)],dims,:))]
    deviation2=[deviation;-deviation[2:end]]
    deviation2=deviation2/first(findmax(abs.(vcat(deviation2...))))

    # create periodic solution with with double period
    new_profile = [psol1.profile; psol1.profile[2:end]] .+ 0.001deviation2
    new_mesh = [psol1.mesh; psol1.mesh[2:end] .+ 1.0] ./ 2
    new_period = 2psol1.period

    psol_guess = psol(
        new_profile,
        psol1.parameters,
        new_mesh,
        new_period,
        psol1.ncol,
        nothing,
        nothing)

    f = (x,psol1) -> vcat([
        psol_res(jet,
            psol([reshape(x[1:end-2],2,:)[:,j] for j in 1:length(psol1.mesh)],
            [psol1.parameters[1:con_par-1]; x[end-1]; psol1.parameters[con_par+1:end]], psol1.mesh, x[end], 
                psol1.ncol, nothing, nothing),psol1,τs)...; 0.0] ...)

    df = (x,psol1) -> defsystem_psol_jac(jet,
            psol([reshape(x[1:end-2],2,:)[:,j] for j in 1:length(psol1.mesh)],
            [psol1.parameters[1:con_par-1]; x[end-1]; psol1.parameters[con_par+1:end]], psol1.mesh, x[end], 
                psol1.ncol, nothing, nothing),psol1,τs)

    x₀ = vec(psol_guess,con_par)
    jac = df(x₀,psol_guess)
    V = [jac[1:end-1,:]; rand(length(x₀))']\[zeros(length(x₀)-1); 1.0]
    x₀new, _, V_new = newton(f, df, x₀, V, psol_guess; tol=1e-8)
    psol_corrected = vec_to_point(x₀new,psol_guess,con_par)

    # continuation with own newton
    psol_branch = (points = psol[],
        tangents = [],
        stepsizes = [],
        f = f,
        df = df,
        parameterbounds = parameterbounds,
        δ=δ,
        δmin=δmin,
        δmax=δmax,
        MaxNumberofSteps = MaxNumberofSteps,
        con_par = con_par,
        NumberOfFails = NumberOfFails)
    push!(psol_branch.points, psol_corrected)
    push!(psol_branch.tangents, V_new)
    push!(psol_branch.stepsizes, 0.0)
    psol_branch
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

function remesh_psol(psol1, new_mesh, new_ncol)
    γ = psol1.profile
    ncol = psol1.ncol
    ts = psol1.mesh 
    γ_new = [[γ[1]]; interpolate.(new_mesh[2:end-1],0.0,Ref(γ),Ref(ts),ncol); [γ[end]]]
    psol(γ_new, psol1.parameters, collect(new_mesh), psol1.period, new_ncol, nothing, nothing)
end

function psolLPC_res(jet,q₁,periodicsolution,psol_ref,τs)
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

function detectSpecialPoints(branch)
    # detect NS bifurcations
    num_unstable = [sum(abs.(p.stability) .> 1.0) for p in branch.points]
    indx_ns = findall(x->x==2,num_unstable[2:end] - num_unstable[1:end-1])

    # detect fold bifurcations
    num_unstable = [sum(abs.(p.stability) .> 1.0 .&& real(p.stability) .> 0.0) for p in branch.points]
    indx_fold = findall(x->abs(x)==1,num_unstable[2:end] - num_unstable[1:end-1])

    # detect pd bifurcations
    num_unstable = [sum(abs.(p.stability) .> 1.0 .&& real(p.stability) .< 0.0) for p in branch.points]
    indx_pd = findall(x->abs(x)==1,num_unstable[2:end] - num_unstable[1:end-1])

    specialpoints = (indx_fold = indx_fold, indx_pd  = indx_pd, indx_ns = indx_ns)
    merge(branch, (specialpoints = specialpoints,))
end
