mutable struct stst
    coords::Vector{Float64}
    parameters::Vector{Float64}
    stability::Union{Vector{ComplexF64},Nothing}
end

function vec(stst1::stst)
    vcat([stst1.coords, stst1.parameters]...)
end

function vec_to_point(::Type{stst},vec,n)
    stst(vec[1:n], [vec[n+1]],nothing)
end


