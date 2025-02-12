function defsystem_psol_jac_no_mod(jet,point,τs)
    M = jet.M

    γ = point.profile
    dims = length(point.profile[1])
    ncol = point.ncol
    ntst = convert(Int,(length(γ)-1)/ncol)

    nodes = first(legendre(ncol))
    T = point.period
    ts = point.mesh
    colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1])/2 .* (nodes .+ 1) for i in 0:ntst-1]...)
    testintervals = ts[1:ncol:end]
    par = point.parameters

    # create extended mesh
    ζ = mod(colpoints[1] - τs[end]/T, 1.0) 
    number_of_meshes = abs(convert(Int,fld(colpoints[1] - τs[end]/T,1.0)))
    interval = searchsortedfirst(testintervals, ζ) .- 1

    mesh_ext = point.mesh[(interval-1)*ncol+1:end-1] .- number_of_meshes
    for i = number_of_meshes:-1:1
        mesh_ext = [mesh_ext; point.mesh[1:end-1] .- (i-1)]
    end
    mesh_ext = [mesh_ext; 1.0]
    γ_ext = [γ[(interval-1)*ncol+1:end-1]; repeat(γ[1:end-1],number_of_meshes,1); [γ[end]]]
    testintervals_ext = mesh_ext[1:ncol:end]

    Jac = zeros(dims*length(colpoints),dims*length(mesh_ext))

    for (i,τ) = enumerate(colpoints)
        ζs = τ .- τs/T
        intervals = searchsortedfirst.(Ref(testintervals_ext), ζs) .- 1

        xs = [mesh_ext[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
        γζ = [ γ_ext[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
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
    end

    Jac
end

# taken from DDE-BifTool
function multipliers(jet, point, τs)
    jac = defsystem_psol_jac_no_mod(jet,point,τs)
    s1,s2=size(jac)
    n_ext=s2-s1
    n2_ext=max(n_ext-s1,0)
    B=blockdiag(sparse(1.0I, n2_ext, n2_ext),sparse(jac[1:s1,n_ext+1:s2]))
    Atop=[spzeros(n2_ext,s1) sparse(1.0I, n2_ext, n2_ext)]
    A=[Atop[:,1:n_ext];-jac[1:s1,1:n_ext]]
    Ptau=[spzeros(n_ext,max(0,2*s1-s2)) sparse(1.0I, n_ext, n_ext)]

    # F = lu(B)
    # Qm = I(length(F.q))[:,F.q]
    # Qm[end-s1+1:end,:]*(F.U\(F.L\((F.Rs .* (A*sparse(1.0I, n_ext, n_ext)))[F.p, :])))
    # (B\Matrix(A*sparse(1.0I, n_ext, n_ext)))[end-s1+1:end,:]

    vals = first(eigen(Ptau*(B\Matrix(A*sparse(1.0I, n_ext, n_ext))))) # fix me: there shouln't be Matrix here!
    # remove trivial multiplier
    deleteat!(vals,findmin(abs.(vals .- 1.0))[2])
    @set point.stability = vals 
end

function multipliers!(jet, branch::Branch, τs)
    for (i,p) in enumerate(branch.points)
        branch.points[i] = multipliers(jet,p,τs)
    end
end
