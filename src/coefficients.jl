function nmfm_coefficient(jet, periodicsolution::PsolLPC, τs, ap)

  T = periodicsolution.period
  par = periodicsolution.parameters[ap]
  ts = T * vec(periodicsolution.mesh)
  γ = periodicsolution.profile

  M = jet.M
  γ = periodicsolution.profile
  dims = length(γ[1])
  ncol = periodicsolution.ncol
  ntst = convert(Int, (length(γ) - 1) / ncol)

  nodes = first(legendre(ncol))
  colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1]) / 2 .* (nodes .+ 1) for i in 0:ntst-1]...)
  testintervals = ts[1:ncol:end]

  γτs = [hcat(interpolateT.(ζ, -τs, Ref(γ), T, Ref(ts), ncol)...) for ζ ∈ colpoints[:]]
  ζ(φ) = vcat(φ.(colpoints[:], 0) + 
              sum(M[j+1].(γτs, Ref(par)) .* (τ * φ.(colpoints[:], -τ)) 
                  for (j, τ) ∈ enumerate(τs[2:end]))...)
  Bτ(γ, φ₁, φ₂) = vcat(Dγ.(γ, φ₁, φ₂, Ref(jet), colpoints[:], Ref(τs), Ref(par))...)
  dγ(τ, θ) = d_interpolateT(τ, θ, γ, T, ts, ncol)

  jac = differential_equation_part(M, γ, ts, par, τs, T, colpoints, dims, ntst, ncol, testintervals)

  q1, p1 = borderedInverse(jac, [-ζ(dγ); zeros(dims + 1)], ntst, ncol, dims, ts, normalization=false)
  p1 = p1[1:end-2]

  # normalize \int_^T <p1,p1> = 1 with Newton-Cotes
  # weights = newtonCotesWeights(ncol)
  # res = 0.0
  # for j ∈ 0:ntst-1, i ∈ 0:ncol
  #   res += (ts[(j+1)*ncol+1] - ts[j*ncol+1]) * dot(p1[j*ncol+1+i], p1[j*ncol+1+i]) * weights[i+1]
  # end
  # p1 /= sqrt(res)

  γ_interpolate(τ, θ) = interpolateT(τ, θ, γ, T, ts, ncol)
  φ(τ, θ) = interpolateT(τ, θ, q1, T, ts, ncol) + θ * dγ(τ, θ)
  φ2(τ, θ) = interpolateT(τ, θ, q1, T, ts, ncol) + 0.5 * θ * dγ(τ, θ)
  d2φ(τ, θ) = d_interpolateT(τ, θ, q1, T, ts, ncol) + 0.5 * θ * dd_interpolateT(τ, θ, γ, T, ts, ncol)

  # The following expression should vanishing
  # @show dot(p1, vcat(dγ.(colpoints[:],0) + sum(M[j+1].(γτs,Ref(par)) .* ( τ*dγ.(colpoints[:],-τ)) for (j,τ) ∈ enumerate(τs[2:end]))...))

  dot(p1, Bτ(γ_interpolate, φ, φ) - 2.0 * ζ(d2φ)) / (2.0 * dot(p1, ζ(φ2)))
end

function nmfm_coefficient(jet, periodicsolution::psol_ns, τs, ap)
  M = jet.M

  γ = periodicsolution.profile
  dims = length(γ[1])
  ncol = periodicsolution.ncol
  ntst = convert(Int, (length(γ) - 1) / ncol)
  T = periodicsolution.period
  par = periodicsolution.parameters[ap]
  ω = periodicsolution.omega
  ts = T * vec(periodicsolution.mesh)
  nodes = first(legendre(ncol))
  colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1]) / 2 .* (nodes .+ 1) for i in 0:ntst-1]...)
  testintervals = ts[1:ncol:end]

  γτs = [hcat(interpolateT.(ζ, -τs, Ref(γ), T, Ref(ts), ncol)...) for ζ ∈ colpoints[:]]
  # ζ(φ) = vcat(φ.(colpoints[:], 0) + sum(M[j].(γτs, Ref(par)) .* (τs[j] * φ.(colpoints[:], -τs[j])) for j ∈ 2:length(τs))...)
  ζ(φ) = vcat(φ.(colpoints[:], 0) + sum(M[j+1].(γτs, Ref(par)) .* (τ * φ.(colpoints[:], -τ)) for (j, τ) ∈ enumerate(τs[2:end]))...)

  jac = differential_equation_part_complex(M, γ, ts, par, τs, T, im * ω / T, colpoints, dims, ntst, ncol, testintervals)
  # q1, p1 = borderedInverse(jac, [zeros(dims * (ncol * ntst + 1)); 1.0], ntst, ncol, dims, ts, normalization=false)

  _, s, V = svd(jac[1:end-1, :]')
  indxmin = last(findmin(s))
  p1 = V[:, indxmin]

  _, s, V = svd(jac[1:end-1, :])
  indxmin = last(findmin(s))
  q1 = V[:, indxmin]

  # q1 = q1[1:end-1]
  q1 = [vec(reshape(q1[dims*(i-1)+1:dims*i], dims, 1)) for i ∈ eachindex(ts)]
  p1 = p1[1:end-dims]

  # normalize \int_^T <q1,q1> = 1 with Newton-Cotes
  weights = newtonCotesWeights(ncol)
  res = 0.0
  for j ∈ 0:ntst-1, i ∈ 0:ncol
    res += (ts[(j+1)*ncol+1] - ts[j*ncol+1]) * dot(q1[j*ncol+1+i], q1[j*ncol+1+i]) * weights[i+1]
  end
  q1 /= sqrt(res)

  γ_interpolate(τ, θ) = interpolateT(τ, θ, γ, T, ts, ncol)
  φ(τ, θ) = interpolateT(τ, θ, q1, T, ts, ncol) * exp(im * ω * θ / T)
  dφ(τ, θ) = d_interpolateT(τ, θ, q1, T, ts, ncol) * exp(im * ω * θ / T)

  # check if q1 is an eigenfunction
  # ζ0(φ) = vcat(sum(M[j].(γτs, Ref(par)) .* (φ.(colpoints[:], -τ)) for (j, τ) ∈ enumerate(τs))...)
  # norm(vcat(dφ.(colpoints, 0.0)...) - ζ0(φ) + im * ω * vcat(φ.(colpoints, 0.0)...) / T)

  # define multi-linear forms
  Bτ(φ₁, φ₂) = vcat(Dγ.(Ref(γ_interpolate), φ₁, φ₂, Ref(jet), colpoints[:], Ref(τs), Ref(par))...)
  Cτ(φ₁, φ₂, φ₃) = vcat(Dγ.(Ref(γ_interpolate), φ₁, φ₂, φ₃, Ref(jet), colpoints[:], Ref(τs), Ref(par))...)

  rhs = Bτ(φ, φ)
  jac = differential_equation_part_complex(M, γ, ts, par, τs, T, 2im * ω / T, colpoints, dims, ntst, ncol, testintervals)
  # sort(eigen(jac[1:end-1,:]).values, by=x->norm(x))[1:3]
  u20 = jac[1:end-1, :] \ [rhs; zeros(dims)]
  u20 = [vec(reshape(u20[dims*(i-1)+1:dims*i], dims, 1)) for i ∈ eachindex(ts)]
  # norm(jac[1:end-1,:]*vcat(u20...) - rhs)
  h20(τ, θ) = interpolateT(τ, θ, u20, T, ts, ncol) * exp(2im * ω * θ / T)

  #dh20(τ, θ) = d_interpolateT(τ, θ, u20, T, ts, ncol) * exp(2im * ω * θ / T)
  #norm(vcat(dh20.(colpoints, 0.0)...) - ζ0(h20) + 2 * im * ω * vcat(h20.(colpoints, 0.0)...) / T - Bτ(φ, φ))

  jac = differential_equation_part(M, γ, ts, par, τs, T, colpoints, dims, ntst, ncol, testintervals)
  r = leftnullvector(jac, ntst, ncol, dims)
  # _, s, V = svd(jac[1:end-1, :]')
  # indxmin = last(findmin(s))
  # r = V[:, indxmin][1:end-2]

  dγ(τ, θ) = d_interpolateT(τ, θ, γ, T, ts, ncol)

  a = dot(r, Bτ(φ, conj ∘ φ)) / dot(r, ζ(dγ))
  rhs = real(Bτ(φ, conj ∘ φ) - a * ζ(dγ))
  u11, _ = borderedInverse(jac, [rhs; zeros(dims + 1)], ntst, ncol, dims, ts, normalization=false)
  h11(τ, θ) = interpolateT(τ, θ, u11, T, ts, ncol) + θ * a * d_interpolateT(τ, θ, γ, T, ts, ncol)

  #dh11(τ, θ) = d_interpolateT(τ, θ, u11, T, ts, ncol)
  #norm(vcat(dh11.(colpoints, 0.0)...) - ζ0(h11) - real(Bτ(φ, conj ∘ φ) - a * vcat(dγ.(colpoints, 0.0)...)))

  d = 1 / 2 * dot(p1, Cτ(φ, φ, conj ∘ φ) + 2Bτ(φ, h11) + Bτ(h20, conj ∘ φ) - 2 * a * ζ(dφ)) / dot(p1, ζ(φ))

  d
end

function differential_equation_part(M, γ, ts, par, τs, T, colpoints, dims, ntst, ncol, testintervals, sign1=-, sign2=-, sign3=-, advanced=false)
  Jac = zeros(dims * (ntst * ncol + 1) + 1, dims * (ntst * ncol + 1))

  for (i, τ) = enumerate(colpoints)
    ζs = mod1.(sign2.(τ, τs), T)
    intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1
    if sign3 == +
      if advanced
        signs = -sign.(τ .+ τs[2:end] .- T)
      else
        signs = sign.(τ .- τs[2:end])
      end
    else
      signs = ones(length(τs[2:end]))
    end

    if length(unique(intervals)) < length(τs)
      # println("Delay larger than period.")
      # break
    end

    xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    γs = [γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    γζ = [L(x, y, ζ, ncol) for (ζ, x, y) in zip(ζs, xs, γs)]

    if advanced
      M1 = M[1](hcat(γζ...), par)'
    else
      M1 = M[1](hcat(γζ...), par)
    end
    for k = 0:ncol
      Jac[(i-1)*dims+1:i*dims, ncol*dims*(intervals[1]-1).+range(k * dims + 1, (k + 1) * dims)] += dlj(xs[1], ζs[1], k, ncol) * diagm(ones(dims)) + sign1(M1 * l0(xs[1], ζs[1], k, ncol) * diagm(ones(dims)))
    end
    # delay terms
    for j ∈ 2:length(τs)
      if advanced
        γζ = γζs(ζs[j], τs, ncol, T, testintervals, ts, γ)
        Mj = M[j](hcat(γζ...), par)'
      else
        Mj = M[j](hcat(γζ...), par)
      end
      for k = 0:ncol
        Jac[(i-1)*dims+1:i*dims, ncol*dims*(intervals[j]-1).+range(k * dims + 1, (k + 1) * dims)] += signs[j-1] * sign1(Mj * l0(xs[j], ζs[j], k, ncol) * diagm(ones(dims)))
      end
    end
  end
  Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1), 1:dims] = diagm(ones(dims))
  Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1), end-dims+1:end] = sign3(diagm(ones(dims)))

  Jac
end

function differential_equation_part_complex(M, γ, ts, par, τs, T, σ, colpoints, dims, ntst, ncol, testintervals, sign1=-, sign2=-, sign3=-, advanced=false)
  Jac = zeros(ComplexF64, dims * (ntst * ncol + 1) + 1, dims * (ntst * ncol + 1))

  Id = diagm(ones(dims))
  for (i, τ) = enumerate(colpoints)
    ζs = mod1.(sign2.(τ, τs), T)
    intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1
    if sign3 == +
      if advanced
        signs = -sign.(τ .+ τs[2:end] .- T)
      else
        signs = sign.(τ .- τs[2:end])
      end
    else
      signs = ones(length(τs[2:end]))
    end

    xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    γs = [γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    γζ = [L(x, y, ζ, ncol) for (ζ, x, y) in zip(ζs, xs, γs)]

    if advanced
      M1 = M[1](hcat(γζ...), par)'
    else
      M1 = M[1](hcat(γζ...), par)
    end
    for k = 0:ncol
      Jac[(i-1)*dims+1:i*dims, ncol*dims*(intervals[1]-1).+range(k * dims + 1, (k + 1) * dims)] += dlj(xs[1], ζs[1], k, ncol) * Id +
                                                                                                   sign1(M1 * l0(xs[1], ζs[1], k, ncol) * Id) + σ * Id * l0(xs[1], ζs[1], k, ncol)
    end
    # delay terms
    for j ∈ 2:length(τs)
      if advanced
        γζ = γζs(ζs[j], τs, ncol, T, testintervals, ts, γ)
        Mj = M[j](hcat(γζ...), par)'
      else
        Mj = M[j](hcat(γζ...), par)
      end
      for k = 0:ncol
        Jac[(i-1)*dims+1:i*dims, ncol*dims*(intervals[j]-1).+range(k * dims + 1, (k + 1) * dims)] += signs[j-1] * sign1(Mj * exp(-σ * τs[j]) * l0(xs[j], ζs[j], k, ncol) * Id)
      end
    end
  end
  Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1), 1:dims] = Id
  Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1), end-dims+1:end] = sign3(Id)

  Jac
end

function differential_equation_part_complex_new(M, γ, ts, par, τs, T, σ, colpoints, dims, ntst, ncol, testintervals)
  Jac = zeros(ComplexF64, dims * (ntst * ncol + 1) + 1, dims * (ntst * ncol + 1))

  Id = diagm(ones(dims))
  for (i, τ) = enumerate(colpoints)
    ζs = mod1.(-(τ, τs), T)
    intervals = searchsortedfirst.(Ref(testintervals), ζs) .- 1

    xs = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    γs = [γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    γζ = [L(x, y, ζ, ncol) for (ζ, x, y) in zip(ζs, xs, γs)]

    M1 = M[1](hcat(γζ...), par)
    for k = 0:ncol
      jac_range = (i-1)*dims+1:i*dims, ncol * dims * (intervals[1] - 1) .+ range(k * dims + 1, (k + 1) * dims)
      Jac[jac_range] += dlj(xs[1], ζs[1], k, ncol) * Id - (M1 * l0(xs[1], ζs[1], k, ncol) * Id) + σ * Id * l0(xs[1], ζs[1], k, ncol)
    end

    # delay terms
    for j ∈ 2:length(τs)
      Mj = M[j](hcat(γζ...), par)
      for k = 0:ncol
        jac_range = (i-1)*dims+1:i*dims, ncol * dims * (intervals[j] - 1) .+ range(k * dims + 1, (k + 1) * dims)
        Jac[jac_range] += -(Mj * exp(-σ * τs[j]) * l0(xs[j], ζs[j], k, ncol) * Id)
      end
    end
  end
  Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1), 1:dims] = Id
  Jac[dims*ntst*ncol+1:dims*(ntst*ncol+1), end-dims+1:end] = -Id

  Jac
end

function nmfm_coefficient(jet, periodicsolution::psol_pd, τs, ap)

  M = jet.M

  T = periodicsolution.period
  par = periodicsolution.parameters[ap]
  ts = T * vec(periodicsolution.mesh)
  γ = periodicsolution.profile
  ncol = periodicsolution.ncol

  dims = length(γ[1])
  ntst = convert(Int, (length(γ) - 1) / ncol)
  nodes = first(legendre(ncol))
  colpoints = hcat([ts[i*ncol+1] .+ (ts[(i+1)*ncol+1] - ts[i*ncol+1]) / 2 .* (nodes .+ 1) for i in 0:ntst-1]...)
  testintervals = ts[1:ncol:end]

  γτs = [hcat(interpolateT.(ζ, -τs, Ref(γ), T, Ref(ts), ncol)...) for ζ ∈ colpoints[:]]
  ζ(φ) = vcat(φ.(colpoints[:], 0) + sum(M[j].(γτs, Ref(par)) .* (τs[j] * φ.(colpoints[:], -τs[j])) for j ∈ 2:length(τs))...)
  Bτ(γ, φ₁, φ₂) = vcat(Dγ.(γ, φ₁, φ₂, Ref(jet), colpoints[:], Ref(τs), Ref(par))...)
  Cτ(γ, φ₁, φ₂, φ₃) = vcat(Dγ.(γ, φ₁, φ₂, φ₃, Ref(jet), colpoints[:], Ref(τs), Ref(par))...)
  dγ(τ, θ) = d_interpolateT(τ, θ, γ, T, ts, ncol)

  # eigenfunction
  jac = differential_equation_part(M, γ, ts, par, τs, T, colpoints, dims, ntst, ncol, testintervals, -, -, +)
  q, ptilde = borderedInverse(jac, [zeros(dims * (ntst * ncol + 1)); 1], ntst, ncol, dims, ts, normalization=false)
  ptilde = ptilde[1:end-2]

  # normalize \int_^T <q,q> = 1 with Newton-Cotes
  weights = newtonCotesWeights(ncol)
  res = 0.0
  for j ∈ 0:ntst-1, i ∈ 0:ncol
    res += (ts[(j+1)*ncol+1] - ts[j*ncol+1]) * dot(q[j*ncol+1+i], q[j*ncol+1+i]) * weights[i+1]
  end
  q /= sqrt(res)

  φ(τ, θ) = sign(τ + θ) * (iseven(div(τ + θ, T)) * 2 - 1) * interpolateT(τ, θ, q, T, ts, ncol)
  dφ(τ, θ) = sign(τ + θ) * (iseven(div(τ + θ, T)) * 2 - 1) * d_interpolateT(τ, θ, q, T, ts, ncol)

  jac = differential_equation_part(M, γ, ts, par, τs, T, colpoints, dims, ntst, ncol, testintervals, -, -, -, false)
  _, rtilde = borderedInverse(jac, [zeros(dims * (ntst * ncol + 1)); 1], ntst, ncol, dims, ts, normalization=false)

  γ_interpolate(τ, θ) = interpolateT(τ, θ, γ, T, ts, ncol)
  a = dot(rtilde[1:end-2], Bτ(γ_interpolate, φ, φ)) / dot(rtilde[1:end-2], ζ(dγ)) / 2
  B = h2_RHS(jet, M, par, φ, γ, a, ts, τs, T, colpoints, testintervals, dims, ncol)
  u2, _ = borderedInverse(jac, B, ntst, ncol, dims, ts, normalization=false)
  h2(τ, θ) = interpolateT(τ, θ, u2, T, ts, ncol) + 2 * a * θ * dγ(τ, θ)

  c = dot(ptilde, Cτ(γ_interpolate, φ, φ, φ) + 3 * Bτ(γ_interpolate, φ, h2) - 6 * a * ζ(dφ))
  c /= 6 * dot(ptilde, ζ(φ))
  c

  # J₁τ(γ) = vcat(Dγ.(γ, Ref(jet), colpoints[:], Ref(τs), Ref(par))...)

  # K01 = [-dot(rtilde[1:end-2], J₁τ(γ_interpolate) * [0.0; 1.0]);
  #   dot(rtilde[1:end-2], J₁τ(γ_interpolate) * [1.0; 0.0])]

  # dot(rtilde[1:end-2], J₁τ(γ_interpolate) * K01)
end

function h2_RHS(jet, M, par, φ, γ, a, ts, τs, T, colpoints, testintervals, dims, ncol)
  B = [zeros(dims) for _ ∈ eachindex(colpoints)]
  for i = eachindex(colpoints)
    τ = mod1.(colpoints[i] .- τs, T)
    intervals = searchsortedfirst.(Ref(testintervals), τ) .- 1

    x = [ts[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    y = [γ[1+(interval-1)*ncol:1+interval*ncol] for interval in intervals]
    γτ = [L(x, y, τ, ncol) for (x, y, τ) ∈ zip(x, y, τ)]
    dγ = [dL(x, y, τ, ncol) for (x, y, τ) ∈ zip(x, y, τ)]
    B[i] = -dγ[1]
    for j ∈ 2:length(τs)
      B[i] -= τs[j] * M[j](hcat(γτ...), par) * dγ[j]
    end
    B[i] *= (2.0 * a)
    Φ = hcat(φ.(τ[1], -τs)...)
    B[i] += jet.D2(hcat(γτ...), par, Φ, Φ)
  end
  [vcat(B...); zeros(dims + 1)]
end

function normal_form_coefficients!(jet, branch, τs)
  for (i, p) in enumerate(branch.points)
    branch.points[i] = @set p.nmfm = nmfm_coefficient(jet, p, τs, branch.con_par)
  end
end
