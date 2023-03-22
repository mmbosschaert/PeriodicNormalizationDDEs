function locate_double_hopf(branch)

    num_unstable = zeros(Int,length(branch))
    for (i,p) in enumerate(branch)
        p.stability
        # remove ω from stability field
        ind1 = findfirst(λ -> imag(λ) ≈ p.ω, p.stability)
        ind2 = findfirst(λ -> imag(λ) ≈ -p.ω, p.stability)
        eigenvalues = p.stability[1:end .!= ind1 .&& 1:end .!= ind2]
        # count unstable eigenvalues
        num_unstable[i] = sum(real(eigenvalues) .> 0.0)
    end
    ind_double_hopf = findall(p -> abs(p) == 2, num_unstable[2:end] - num_unstable[1:end-1])
    ind_double_hopf
end
