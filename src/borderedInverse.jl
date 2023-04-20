function borderedInverse(jac,rhs,ntst,ncol,dims,ts;normalization=true)
    # jacex = [jac[1:end-1,:] rand(typeof(first(jac)),dims*(ntst*ncol+1))
    #          rand(typeof(first(jac)),1,dims*(ntst*ncol+1)) 0.0]
    # p1 = jacex'\[zeros(dims*(ntst*ncol+1)); 1.0]

    _, s, V = svd(jac[1:end-1,:]')
    indxmin = last(findmin(s))
    p1 = V[:, indxmin]
    _, s, V = svd(jac[1:end-1,:])
    indxmin = last(findmin(s))
    q1 = V[:, indxmin]

    if normalization
        sol = [jac p1]\rhs
    else
        # q1 = jacex\[zeros(dims*(ntst*ncol+1)); 1.0]
        sol = [jac[1:end-1,:] p1; q1' 0.0]\rhs
        # sol = jacex\rhs
    end
    sol = [vec(reshape(sol[dims*(i-1)+1:dims*i],dims,1)) for i ∈ eachindex(ts)]
    sol, p1
end


function leftnullvector(jac,ntst,ncol,dims)
    _, s, V = svd(jac[1:end-1,:]')
    indxmin = last(findmin(s))
    V[:, indxmin][1:end-2]
end
