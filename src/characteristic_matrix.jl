function characteristic_matrices(model, bt, τs)

    # calculate derivatives at Bogdanov-Takens point
    jac = Symbolics.jacobian
    n = length(bt.coords)
    m = length(τs)

    # construct character matrix
    @variables φ[1:n,1:m] α1 α2
    φ = reshape(Symbolics.get_variables(φ), n, m) # cast into a Matrix of symbols
    pars = [α1 α2]

    subs = Dict(φ .=> repeat(bt.coords,1,m))
    subs = merge(subs,Dict(pars .=> bt.pars'))

    @variables z
    Δ = Matrix(z*I,n,n)
    for j=1:length(τs)
        Mⱼ = substitute(jac(model(φ, bt.pars), vcat(φ[:,j])), subs)
        Δ -= Mⱼ*exp(-z*τs[j])
    end 
    Δ = build_function(Δ, z, expression=Val{false})[1]

    # define derivatives of the characteristic matrix evaluated at 0
    d(i) = expand_derivatives∘Differential(z)^i
    Δ′ = build_function(d(1).(Δ(z)), z, expression=Val{false})[1](0)
    Δ′′ = build_function(d(2).(Δ(z)), z, expression=Val{false})[1](0)
    Δ′′′ = build_function(d(3).(Δ(z)), z, expression=Val{false})[1](0)
    Δ⁴ = build_function(d(4).(Δ(z)), z, expression=Val{false})[1](0)
    Δ⁵ = build_function(d(5).(Δ(z)), z, expression=Val{false})[1](0)

    Δ, Δ′, Δ′′, Δ′′′, Δ⁴, Δ⁵

end
