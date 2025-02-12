# function coefficient_matrices(model, point, τs)
function coefficient_matrices(model, n, τs)

    # calculate derivatives at  point
    jac = Symbolics.jacobian
    # n = length(point.coords)
    m = length(τs)

    # construct character matrix
    @variables φ[1:n,1:m] α1 α2
    φ = reshape(Symbolics.get_variables(φ), n, m) # cast into a Matrix of symbols
    pars = [α1 α2]

    M = []
    for j=1:length(τs)
        # push!(M, (coords,pars) -> build_function(jac(model(φ, pars), vcat(φ[:,j])), φ,pars,expression=Val{false})[1](repeat(coords,m),pars))
        push!(M, jac(model(φ, pars), vcat(φ[:,j])))
    end 

    build_function(M,φ,pars,expression=Val{false})[1]
end
