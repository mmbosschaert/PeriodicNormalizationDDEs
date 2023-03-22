function stability(branch,M,τs)
    for p in branch
        dep = DEP(M(repeat(p.coords,1,length(τs)),p.parameters), vec(τs))
        λ,~ = iar_chebyshev(dep,maxit=100,neigs=13)
        p.stability = λ
    end
end
