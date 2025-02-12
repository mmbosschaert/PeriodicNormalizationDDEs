function l0(x,ζ,j,n)
    Π = 1.0
    @inbounds for m ∈ 0:n 
        if m != j
            Π *= (ζ-x[m+1])/(x[j+1]-x[m+1])
        end
    end
    Π
end

function dlj(x,ζ,j,n)
    S = 0.0;
    @inbounds for i=0:n
        if i != j
            Π = 1.0
            for m ∈ 0:n
                if m != i && m != j
                    Π *= (ζ-x[m+1])/(x[j+1]-x[m+1])
                end
            end
            S += 1/(x[j+1]-x[i+1])*Π
        end
    end
    S
end

function ddlj(x,ζ,j,n)
    S = 0.0;
    for i=0:n
        if i != j
            S2 = 0.0;
            for m=0:n
                if m != i && m != j
                    Π = 1.0
                    for l ∈ 0:n
                        if l != i && l != j && l != m
                            Π *= (ζ-x[l+1])/(x[j+1]-x[l+1])
                        end
                    end
                    S2  += 1/(x[j+1]-x[m+1])*Π
                end
            end
            S  += 1/(x[j+1]-x[i+1])*S2
        end
    end
    S
end

function L(x,y,ζ,n)
    vals = zeros(size(y[1]))
    for j ∈ 0:n
        lj = l0(x,ζ,j,n)
        vals += y[j+1].*lj
    end
    vals
end

function L(x,y::Matrix{Float64},ζ,n)
    vals = zeros(size(y,1))
    for j ∈ 0:n
        lj = l0(x,ζ,j,n)
        vals += y[:,j+1].*lj
    end
    vals
end

function dL(x,y,ζ,n)
    val = zeros(size(y[1]))
    for j ∈ 0:n
        temp = dlj(x,ζ,j,n)
        val += y[j+1].*temp
    end
    val
end

function dL(x,y::Matrix{Float64},ζ,n)
    val = zeros(size(y,1))
    for j ∈ 0:n
        temp = dlj(x,ζ,j,n)
        val += y[:,j+1].*temp
    end
    val
end

function ddL(x,y,ζ,n)
    val = zeros(size(y[1]))
    for j ∈ 0:n
        temp = ddlj(x,ζ,j,n)
        val += y[j+1].*temp
    end
    val
end

function interpolate(τ,θ,y,ts,ncol)
        ζ = mod(τ + θ,1.0)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[1+(interval-1)*ncol:1+interval*ncol] 
        L(x,ys,ζ,ncol)
end

function interpolate(τ,θ,y::Matrix{Float64},ts,ncol)
        ζ = mod(τ + θ,1.0)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[:,1+(interval-1)*ncol:1+interval*ncol] 
        L(x,ys,ζ,ncol)
end

function d_interpolate(τ,θ,y,ts,ncol)
        ζ = mod(τ + θ,1.0)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[1+(interval-1)*ncol:1+interval*ncol] 
        vec(dL(x,ys,ζ,ncol))
end

function d_interpolate(τ,θ,y::Matrix{Float64},ts,ncol)
        ζ = mod(τ + θ,1.0)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[:,1+(interval-1)*ncol:1+interval*ncol] 
        vec(dL(x,ys,ζ,ncol))
end

function interpolate(τ,θ,y,T,ts,ncol)
        ζ = mod1(τ + θ/T,1)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[1+(interval-1)*ncol:1+interval*ncol] 
        L(x,ys,ζ,ncol)
end

function d_interpolate(τ,θ,y,T,ts,ncol)
        ζ = mod1(τ + θ/T,1)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[1+(interval-1)*ncol:1+interval*ncol] 
        vec(dL(x,ys,ζ,ncol))
end

function interpolateT(τ,θ,y,T,ts,ncol)
        ζ = mod1(τ + θ,T)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[1+(interval-1)*ncol:1+interval*ncol] 
        L(x,ys,ζ,ncol)
end

function d_interpolateT(τ,θ,y,T,ts,ncol)
        ζ = mod1(τ + θ,T)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[1+(interval-1)*ncol:1+interval*ncol] 
        vec(dL(x,ys,ζ,ncol))
end

function dd_interpolateT(τ,θ,y,T,ts,ncol)
        ζ = mod1(τ + θ,T)
        testintervals = ts[1:ncol:end]
        interval = searchsortedfirst(testintervals, ζ) .- 1
        x = ts[1+(interval-1)*ncol:1+interval*ncol]
        ys = y[1+(interval-1)*ncol:1+interval*ncol] 
        vec(ddL(x,ys,ζ,ncol))
end

# @register_symbolic searchsortedfirst(testintervals::Vector{Float64}, ζ)

# @register_symbolic   interpolate(τ::Float64,θ::Float64,y,T,ts::Vector{Float64},ncol::Int)
# @register_symbolic d_interpolate(τ::Float64,θ::Float64,y,T,ts::Vector{Float64},ncol::Int)

#Symbolics.derivative(::typeof(interpolate), args::NTuple{1,Any}, ::Val{1}) = cos(args[1])
