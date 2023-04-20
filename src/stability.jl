function stability(M,points,τs; neigs = 13, maxit = 100)
    for (i,p) in enumerate(points)
        dep = DEP(M(repeat(p.coords,1,length(τs)),p.parameters), vec(τs))
        λ,_ = iar_chebyshev(dep,maxit=maxit,neigs=neigs)
        points[i] = @set p.stability = λ
    end
end

function stability(jet::NamedTuple,p,τs; neigs = 13, maxit = 100)
    dep = DEP(jet.Ms(repeat(p.coords,1,length(τs)),p.parameters), vec(τs))
    λ,_ = iar_chebyshev(dep,maxit=maxit,neigs=neigs)
    @set p.stability = λ
end
