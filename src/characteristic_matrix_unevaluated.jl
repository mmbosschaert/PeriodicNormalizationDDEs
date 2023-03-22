# this function is not necessary anymore, see getJet
function characteristic_matrices_unevaluated(model, n, τs)

    # calculate derivatives at Bogdanov-Takens point
    jac = Symbolics.jacobian
    m = length(τs)

    # construct character matrix
    @variables φ[1:n,1:m] α1 α2
    φ = reshape(Symbolics.get_variables(φ), n, m) # cast into a Matrix of symbols
    pars = [α1 α2]

    @variables z
    Δ = Matrix(z*I,n,n)
    for j in eachindex(τs)
        Mⱼ = jac(model(φ, pars), vcat(φ[:,j]))
        Δ -= Mⱼ*exp(-z*τs[j])
    end 
    Δ = build_function(Δ, φ, pars, z, expression=Val{false})[1]

end

function characteristic_matrices_unevaluated_re_im(model, n, τs)

    # calculate derivatives at Bogdanov-Takens point
    jac = Symbolics.jacobian
    m = length(τs)

    # construct character matrix
    @variables φ[1:n,1:m] α1 α2
    φ = reshape(Symbolics.get_variables(φ), n, m) # cast into a Matrix of symbols
    pars = [α1 α2]

    @variables z
    Δre = Matrix(0.0*I,n,n)
    Δim = Matrix(z*I,n,n)
    for j in eachindex(τs)
        Mⱼ = jac(model(φ, pars), vcat(φ[:,j]))
        Δre -= Mⱼ*cos(z*τs[j])
        Δim += Mⱼ*sin(z*τs[j])
    end 
    Δre = build_function(Δre, φ, pars, z, expression=Val{false})[1]
    Δim = build_function(Δim, φ, pars, z, expression=Val{false})[1]

    Δre, Δim
end
